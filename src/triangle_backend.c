// Keep the baseline Triangle source and arithmetic unchanged. The C boundary
// converts its exit() paths into a failed call and frees partial allocations.
// Calls are serialized by triangle.cpp; no longjmp crosses a C++ stack frame.
#include <setjmp.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

struct Allocation {
  void* pointer;
  struct Allocation* next;
};
static struct Allocation* allocations;
static jmp_buf triangle_failure;

static void TriangleExit(int status) {
  longjmp(triangle_failure, status + 1);
}
static void* TriangleAllocate(size_t size) {
  void* pointer = malloc(size);
  struct Allocation* record;
  if (!pointer) {
    TriangleExit(1);
  }
  record = (struct Allocation*)malloc(sizeof(*record));
  if (!record) {
    free(pointer);
    TriangleExit(1);
  }
  record->pointer = pointer;
  record->next = allocations;
  allocations = record;
  return pointer;
}
static void TriangleFree(void* pointer) {
  struct Allocation** link = &allocations;
  while (*link && (*link)->pointer != pointer) {
    link = &(*link)->next;
  }
  if (*link) {
    struct Allocation* record = *link;
    *link = record->next;
    free(record);
  }
  free(pointer);
}
#define exit TriangleExit
#define malloc TriangleAllocate
#define free TriangleFree
#define triunsuitable TriangleUnsuitable
#define lucid_auto_stats triangle_auto_stats
#include "../third_party/triangle/triangle.c"
#undef triunsuitable
#undef lucid_auto_stats
#undef exit
#undef malloc
#undef free

int TriangleExecute(char* switches, struct triangulateio* input,
                    struct triangulateio* output) {
  int status = setjmp(triangle_failure);
  if (!status) {
    triangulate(switches, input, output, NULL);
  }
  while (allocations) {
    struct Allocation* record = allocations;
    allocations = record->next;
    if (status) {
      free(record->pointer);
    }
    free(record);
  }
  if (status) {
    memset(output, 0, sizeof(*output));
  }
  return status;
}
