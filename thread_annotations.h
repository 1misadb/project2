#ifndef THREAD_ANNOTATIONS_H
#define THREAD_ANNOTATIONS_H
#ifdef __clang__
#define THREAD_GUARDED_BY(x) __attribute__((guarded_by(x)))
#else
#define THREAD_GUARDED_BY(x)
#endif
#endif
