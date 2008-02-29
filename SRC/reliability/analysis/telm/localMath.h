// Define 32 bit signed and unsigned integers.
// Change these definitions, if necessary, on 64 bit computers

#ifndef _WIN32

#define int32 int
#define uint32 unsigned int

uint32 _lrotl (uint32 x, int r) {
  return (x << r) | (x >> (sizeof(x)*8-r));
}

#endif
