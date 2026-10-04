#include <array>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "tessellation/tessellation.h"

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>

#include <csignal>
#endif

namespace {
namespace fs = std::filesystem;

void Require(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

std::string ReadFile(const fs::path& path) {
  std::ifstream input(path, std::ios::binary);
  Require(bool(input), "Cannot read " + path.u8string());
  return {std::istreambuf_iterator<char>(input),
          std::istreambuf_iterator<char>()};
}

struct ProcessResult {
  int status;
  std::string output;
};

#ifdef _WIN32
struct Handle {
  HANDLE value;
  explicit Handle(HANDLE handle) : value(handle) {}
  ~Handle() {
    if (value && value != INVALID_HANDLE_VALUE) {
      CloseHandle(value);
    }
  }
  Handle(const Handle&) = delete;
  Handle& operator=(const Handle&) = delete;
};

std::wstring Wide(const std::string& text) {
  if (text.empty()) {
    return {};
  }
  const int size =
      MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, text.data(),
                          static_cast<int>(text.size()), nullptr, 0);
  Require(size > 0, "Invalid UTF-8 process argument");
  std::wstring result(size, L'\0');
  MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, text.data(),
                      static_cast<int>(text.size()), result.data(), size);
  return result;
}

std::wstring Quote(const std::wstring& argument) {
  std::wstring result = L"\"";
  std::size_t backslashes = 0;
  for (const wchar_t character : argument) {
    if (character == L'\\') {
      ++backslashes;
      continue;
    }
    result.append(character == L'"' ? 2 * backslashes + 1 : backslashes, L'\\');
    result.push_back(character);
    backslashes = 0;
  }
  result.append(2 * backslashes, L'\\');
  result.push_back(L'"');
  return result;
}

std::vector<wchar_t> StandaloneEnvironment(const fs::path& directory) {
  const DWORD size = GetEnvironmentVariableW(L"SystemRoot", nullptr, 0);
  Require(size > 0, "SystemRoot is missing");
  std::wstring system_root(size, L'\0');
  system_root.resize(
      GetEnvironmentVariableW(L"SystemRoot", system_root.data(), size));
  const std::map<std::wstring, std::wstring> variables{
      {L"PATH", system_root + L"\\System32"},
      {L"SystemRoot", system_root},
      {L"TEMP", directory.wstring()},
      {L"TMP", directory.wstring()},
      {L"WINDIR", system_root}};
  std::vector<wchar_t> environment;
  for (const auto& variable : variables) {
    const auto entry = variable.first + L"=" + variable.second;
    environment.insert(environment.end(), entry.begin(), entry.end());
    environment.push_back(L'\0');
  }
  environment.push_back(L'\0');
  return environment;
}
#endif

ProcessResult RunProcess(const fs::path& executable,
                         const std::vector<std::string>& arguments,
                         const fs::path& directory) {
  const auto log = directory / "process.log";
  int status;
#ifdef _WIN32
  {
    SECURITY_ATTRIBUTES security{sizeof(SECURITY_ATTRIBUTES), nullptr, TRUE};
    Handle output(CreateFileW(log.c_str(), GENERIC_WRITE, FILE_SHARE_READ,
                              &security, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL,
                              nullptr));
    Require(output.value != INVALID_HANDLE_VALUE,
            "Cannot create the process log");
    std::wstring command = Quote(executable.wstring());
    for (const auto& argument : arguments) {
      command += L" " + Quote(Wide(argument));
    }
    STARTUPINFOW startup{};
    startup.cb = sizeof(startup);
    startup.dwFlags = STARTF_USESTDHANDLES;
    startup.hStdOutput = startup.hStdError = output.value;
    PROCESS_INFORMATION process{};
    std::vector<wchar_t> environment;
#if !defined(_DEBUG)
    // Release uses the installed VC runtime and adjacent iko DLL, with no
    // OCCT or developer-tool paths. Debug needs the nonredistributable CRT.
    environment = StandaloneEnvironment(directory);
#endif
    if (!CreateProcessW(executable.c_str(), command.data(), nullptr, nullptr,
                        TRUE, CREATE_NO_WINDOW | CREATE_UNICODE_ENVIRONMENT,
                        environment.empty() ? nullptr : environment.data(),
                        directory.c_str(), &startup, &process)) {
      throw std::runtime_error("Cannot start the CLI: " +
                               std::to_string(GetLastError()));
    }
    Handle process_handle(process.hProcess);
    Handle thread_handle(process.hThread);
    if (WaitForSingleObject(process_handle.value, 90000) != WAIT_OBJECT_0) {
      TerminateProcess(process_handle.value, 1);
      WaitForSingleObject(process_handle.value, 5000);
      throw std::runtime_error("The CLI did not finish within 90 seconds");
    }
    DWORD exit_code;
    Require(GetExitCodeProcess(process_handle.value, &exit_code) != 0,
            "Cannot read the CLI exit code");
    status = static_cast<int>(exit_code);
  }
#else
  const pid_t child = fork();
  Require(child >= 0, "Cannot start the CLI");
  if (child == 0) {
    const int output = open(log.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0600);
    if (output < 0 || chdir(directory.c_str()) != 0) {
      _exit(127);
    }
    dup2(output, STDOUT_FILENO);
    dup2(output, STDERR_FILENO);
    close(output);
    std::vector<std::string> values{executable.string()};
    values.insert(values.end(), arguments.begin(), arguments.end());
    std::vector<char*> pointers;
    for (auto& value : values) {
      pointers.push_back(value.data());
    }
    pointers.push_back(nullptr);
    execv(executable.c_str(), pointers.data());
    _exit(127);
  }
  const auto deadline =
      std::chrono::steady_clock::now() + std::chrono::seconds(90);
  int wait_status;
  while (waitpid(child, &wait_status, WNOHANG) == 0) {
    if (std::chrono::steady_clock::now() >= deadline) {
      kill(child, SIGKILL);
      waitpid(child, &wait_status, 0);
      throw std::runtime_error("The CLI did not finish within 90 seconds");
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
  }
  status = WIFEXITED(wait_status) ? WEXITSTATUS(wait_status)
                                  : 128 + WTERMSIG(wait_status);
#endif
  return {status, ReadFile(log)};
}

std::uint32_t Word(const std::string& bytes, std::size_t offset) {
  std::uint32_t value = 0;
  for (std::size_t byte = 0; byte < 4; ++byte) {
    value |= std::uint32_t(static_cast<unsigned char>(bytes.at(offset + byte)))
             << (8 * byte);
  }
  return value;
}

float Scalar(const std::string& bytes, std::size_t offset) {
  const auto bits = Word(bytes, offset);
  float value;
  static_assert(sizeof(value) == sizeof(bits));
  std::memcpy(&value, &bits, sizeof(value));
  return value;
}

void CheckObj(const fs::path& path, const iko::Mesh& mesh) {
  std::ifstream input(path);
  Require(bool(input), "OBJ output is missing");
  std::vector<iko::Vec3> positions;
  std::vector<iko::Vec2> uv;
  std::vector<iko::Vec3> normals;
  std::vector<iko::Triangle> triangles;
  std::string line;
  while (std::getline(input, line)) {
    std::istringstream fields(line);
    std::string kind;
    fields >> kind;
    if (kind == "v" || kind == "vn") {
      iko::Vec3 value;
      fields >> value[0] >> value[1] >> value[2];
      Require(bool(fields), "Malformed OBJ vertex record");
      (kind == "v" ? positions : normals).push_back(value);
    } else if (kind == "vt") {
      iko::Vec2 value;
      fields >> value[0] >> value[1];
      Require(bool(fields), "Malformed OBJ UV record");
      uv.push_back(value);
    } else if (kind == "f") {
      iko::Triangle triangle;
      for (auto& vertex : triangle) {
        std::string token;
        fields >> token;
        Require(bool(fields), "Malformed OBJ face record");
        vertex = static_cast<std::uint32_t>(
            std::stoul(token.substr(0, token.find('/'))) - 1);
      }
      triangles.push_back(triangle);
    }
  }
  Require(positions == mesh.xyz && uv == mesh.uv && normals == mesh.normals,
          "OBJ vertex arrays differ from the SDK result");
  Require(triangles == mesh.triangles,
          "OBJ directed triangles differ from the SDK result");
}

void CheckStl(const fs::path& path, const iko::Mesh& mesh) {
  const auto bytes = ReadFile(path);
  Require(bytes.size() >= 84 && Word(bytes, 80) == mesh.triangles.size() &&
              bytes.size() == 84 + 50 * mesh.triangles.size(),
          "STL record count or file length is incorrect");
  for (std::size_t face = 0; face < mesh.triangles.size(); ++face) {
    const auto offset = 84 + 50 * face;
    for (std::size_t corner = 0; corner < 3; ++corner) {
      const auto& vertex = mesh.xyz.at(mesh.triangles[face][corner]);
      for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
        Require(Scalar(bytes, offset + 12 + 4 * (3 * corner + coordinate)) ==
                    static_cast<float>(vertex[coordinate]),
                "STL directed facet coordinates differ from the SDK result");
      }
    }
    Require(bytes[offset + 48] == 0 && bytes[offset + 49] == 0,
            "STL attribute count is nonzero");
  }
}

int Run(int argc, char** argv) {
  Require(argc == 5,
          "CLI tests require the executable, data/output directories and DLL");
  const auto executable = fs::absolute(fs::u8path(argv[1]));
  const auto data = fs::absolute(fs::u8path(argv[2]));
  const auto run_id =
      std::chrono::steady_clock::now().time_since_epoch().count();
  const auto out =
      fs::absolute(fs::u8path(argv[3])) / ("cli-" + std::to_string(run_id));
  fs::create_directories(out);
  const auto binary = out / executable.filename();
  fs::copy_file(executable, binary);
  const auto sdk_dll = fs::absolute(fs::u8path(argv[4]));
  fs::copy_file(sdk_dll, out / sdk_dll.filename());
  fs::copy_file(data / "two_planes.igs", out / fs::u8path(u8"模型.igs"));
  fs::copy_file(data / "prismatic_brep.igs", out / "brep.igs");
  std::ofstream(out / "invalid.igs") << "invalid";
  std::size_t checks = 0;
  const auto run = [&](const std::string& name,
                       const std::vector<std::string>& arguments,
                       int expected_status) {
    const auto child = RunProcess(binary, arguments, out);
    Require(child.status == expected_status,
            name + ": exit " + std::to_string(child.status) + ", expected " +
                std::to_string(expected_status) + "\n" + child.output);
    ++checks;
    return child.output;
  };

  Require(run("help", {"--help"}, 0).find("Exit codes:") != std::string::npos,
          "CLI help must describe the exit codes");
  const std::vector<std::pair<std::string, std::vector<std::string>>> invalid{
      {"missing_geometry", {}},
      {"missing_file", {"missing.igs"}},
      {"invalid_file", {"invalid.igs"}},
      {"unknown_option", {u8"模型.igs", "--unknown", "1"}},
      {"missing_option_value", {u8"模型.igs", "--max-edge"}},
      {"invalid_numeric", {u8"模型.igs", "--chord", "nan"}},
      {"zero_chord", {u8"模型.igs", "--chord", "0"}},
      {"min_above_max", {u8"模型.igs", "--min-edge", "1", "--max-edge", ".2"}},
      {"negative_limit", {u8"模型.igs", "--steiner-limit", "-1"}},
      {"fractional_limit", {u8"模型.igs", "--steiner-limit", ".5"}},
      {"missing_surface", {u8"模型.igs", "--surface", "999999"}},
      {"removed_profile", {u8"模型.igs", "--import-profile", "source"}},
      {"removed_assembly", {u8"模型.igs", "--assembly", "independent_blocks"}},
      {"unknown_output", {u8"模型.igs", "--out", "mesh.txt"}},
      {"two_inputs", {u8"模型.igs", u8"模型.igs"}},
      {"invalid_list_parameters", {u8"模型.igs", "--list", "--chord", "0"}}};
  for (const auto& item : invalid) {
    Require(run(item.first, item.second, 2).find("input_problem:") !=
                std::string::npos,
            item.first + ": missing input error diagnostic");
  }
  const auto listed = run("list", {u8"模型.igs", "--list"}, 0);
  std::smatch selected;
  Require(
      std::regex_search(listed, selected, std::regex("(?:^|\n)DE([0-9]+) ")),
      "CLI listing did not return a surface entry");
  run("select_surface", {u8"模型.igs", "--surface", selected[1].str()}, 0);

  const auto expected = [&](const fs::path& input,
                            const iko::Options& options = iko::Options{}) {
    auto result = iko::Tessellate(input, options);
    Require(result.Succeeded() && result.complete && result.whole_model,
            "SDK comparison result is incomplete: " + result.message);
    return result;
  };
  const auto result = expected(out / fs::u8path(u8"模型.igs"));
  Require(result.model.required_faces == result.faces.size() &&
              result.faces.size() == 2,
          "The default CLI case must cover the complete multi-face model");
  const auto summary = run("default_whole_model_unicode", {u8"模型.igs"}, 0);
  Require(summary.find(std::to_string(result.mesh.xyz.size()) + " vertices, " +
                       std::to_string(result.mesh.triangles.size()) +
                       " triangles") != std::string::npos &&
              summary.find("unused vertices") != std::string::npos &&
              summary.find("degenerate facets") != std::string::npos,
          "CLI summary differs from the SDK result");
  for (const auto& entry : fs::directory_iterator(out)) {
    Require(entry.path().extension() != ".obj" &&
                entry.path().extension() != ".stl",
            "The CLI exported a mesh without --out");
  }
  run("obj_export", {u8"模型.igs", "--out", u8"网格.obj"}, 0);
  run("stl_export", {u8"模型.igs", "--out", u8"网格.stl"}, 0);
  CheckObj(out / fs::u8path(u8"网格.obj"), result.mesh);
  CheckStl(out / fs::u8path(u8"网格.stl"), result.mesh);
  run("export_io_failure", {u8"模型.igs", "--out", "absent/mesh.obj"}, 5);
  iko::Options options;
  options.surface.max_edge = 2.5;
  const auto edge = expected(out / fs::u8path(u8"模型.igs"), options);
  run("separate_curve_surface_parameters",
      {u8"模型.igs", "--surface-max", "2.5", "--out", "edge.obj"}, 0);
  CheckObj(out / "edge.obj", edge.mesh);
  const auto brep = expected(out / "brep.igs");
  Require(brep.model.required_faces == brep.faces.size() &&
              brep.faces.size() == 10 && brep.model.solids == 1,
          "The B-rep CLI case must cover all ten solid faces");
  run("brep_complete", {"brep.igs", "--out", "brep.obj"}, 0);
  run("brep_stl", {"brep.igs", "--out", "brep.stl"}, 0);
  CheckObj(out / "brep.obj", brep.mesh);
  CheckStl(out / "brep.stl", brep.mesh);
  Require(checks == 26, "All 26 existing CLI checks must run");
  std::cout << checks << " CLI and export checks passed\n";
  return 0;
}
}  // namespace

int main(int argc, char** argv) {
#ifdef _WIN32
  SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX | SEM_NOOPENFILEERRORBOX);
#endif
  try {
    return Run(argc, argv);
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
