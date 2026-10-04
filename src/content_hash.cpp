#include "tessellation/private/internal.h"

#include <stdexcept>
#include <string>
#include <vector>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>  // IWYU pragma: keep

#include <bcrypt.h>
#include <minwindef.h>
#endif

namespace iko::internal {
std::string ContentSha256(const std::string& bytes) {
#ifdef _WIN32
  BCRYPT_ALG_HANDLE alg = nullptr;
  BCRYPT_HASH_HANDLE hash = nullptr;
  if (BCryptOpenAlgorithmProvider(&alg, BCRYPT_SHA256_ALGORITHM, nullptr, 0) <
      0) {
    throw std::runtime_error("SHA256 initialization failed");
  }
  DWORD count = 0;
  DWORD size = 0;
  BCryptGetProperty(alg, BCRYPT_OBJECT_LENGTH, reinterpret_cast<PUCHAR>(&size),
                    sizeof(size), &count, 0);
  std::vector<unsigned char> object(size);
  std::vector<unsigned char> digest(32);
  auto status =
      BCryptCreateHash(alg, &hash, object.data(), size, nullptr, 0, 0);
  if (status >= 0) {
    status = BCryptHashData(
        hash, reinterpret_cast<PUCHAR>(const_cast<char*>(bytes.data())),
        static_cast<ULONG>(bytes.size()), 0);
  }
  if (status >= 0) {
    status = BCryptFinishHash(hash, digest.data(),
                              static_cast<ULONG>(digest.size()), 0);
  }
  if (hash) {
    BCryptDestroyHash(hash);
  }
  BCryptCloseAlgorithmProvider(alg, 0);
  if (status < 0) {
    throw std::runtime_error("SHA256 failed");
  }
  static const char kHexDigits[] = "0123456789abcdef";
  std::string result;
  for (unsigned char digest_byte : digest) {
    result += kHexDigits[digest_byte >> 4];
    result += kHexDigits[digest_byte & 15];
  }
  return result;
#else
  throw Unsupported(
      "This release targets the baseline Windows x64 environment");
#endif
}
}  // namespace iko::internal
