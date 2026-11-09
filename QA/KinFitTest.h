//****************************************************************************
//*                    This file is part of KinFit.                          *
//*                                                                          *
//*             KinFit is distributed under the terms of the                 *
//*              GNU General Public License (GPL) version 3,                 *
//*                 copied verbatim in the file "LICENSE".                   *
//*                                                                          *
//*                             Copyright 2024                               *
//*             GSI Helmholtzzentrum für Schwerionenforschung                *
//*          This software is distributed under the terms of the             *
//*          GNU General Public Licence version 3 (GPL Version 3)            *
//*                                                                          *
//*     The copyright holders are listed in the file "COPYRIGHTHOLDERS".     *
//*               The authors are listed in the file "AUTHORS".              *
//****************************************************************************

/**
 * KinFitTest.h
 *
 * Minimal header-only test harness for the KinFit QA suite (no external
 * framework). A test executable defines cases with KF_TEST and ends with
 * KF_TEST_MAIN(). Usage:
 *
 *   ./test_foo            run all cases
 *   ./test_foo <case>     run one case (this is what CTest does)
 *   ./test_foo --list     list case names
 *
 * Checks print file:line and continue; a case fails if any check failed.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <functional>
#include <sstream>
#include <string>
#include <vector>

namespace kftest {

struct TestCase {
  std::string name;
  std::function<void()> fn;
};

inline std::vector<TestCase> &registry() {
  static std::vector<TestCase> cases;
  return cases;
}

inline int &failedChecks() {
  static int n = 0;
  return n;
}

struct Registrar {
  Registrar(const char *name, std::function<void()> fn) {
    registry().push_back({name, std::move(fn)});
  }
};

inline void reportFailure(const char *file, int line, const std::string &msg) {
  ++failedChecks();
  std::printf("  FAILED %s:%d: %s\n", file, line, msg.c_str());
}

/** Print an informational line (shown in `ctest --output-on-failure`). */
template <typename... Args> void info(const char *fmt, Args... args) {
  std::printf("  ");
  std::printf(fmt, args...);
  std::printf("\n");
}

inline int runCase(const TestCase &tc) {
  failedChecks() = 0;
  std::printf("[ RUN  ] %s\n", tc.name.c_str());
  try {
    tc.fn();
  } catch (const std::exception &e) {
    reportFailure(__FILE__, __LINE__,
                  std::string("uncaught exception: ") + e.what());
  }
  const bool ok = failedChecks() == 0;
  std::printf("[ %s ] %s\n", ok ? " OK " : "FAIL", tc.name.c_str());
  return ok ? 0 : 1;
}

inline int runTests(int argc, char **argv) {
  if (argc > 1 && std::strcmp(argv[1], "--list") == 0) {
    for (const auto &tc : registry())
      std::printf("%s\n", tc.name.c_str());
    return 0;
  }
  int nFailed = 0, nRun = 0;
  for (const auto &tc : registry()) {
    if (argc > 1 && tc.name != argv[1])
      continue;
    ++nRun;
    nFailed += runCase(tc);
  }
  if (nRun == 0) {
    std::printf("no test case named '%s'\n", argc > 1 ? argv[1] : "");
    return 2;
  }
  return nFailed == 0 ? 0 : 1;
}

} // namespace kftest

#define KF_TEST(name)                                                          \
  static void name();                                                          \
  static kftest::Registrar kf_registrar_##name(#name, name);                   \
  static void name()

#if defined(__has_feature)
#if __has_feature(address_sanitizer)
#define KF_ASAN 1
#endif
#endif
#if defined(__SANITIZE_ADDRESS__)
#define KF_ASAN 1
#endif

// Under AddressSanitizer, ROOT's TBB teardown crashes at exit (macOS); leave
// via _Exit so that only the test result counts.
#ifdef KF_ASAN
#define KF_TEST_MAIN()                                                         \
  int main(int argc, char **argv) {                                            \
    const int rc = kftest::runTests(argc, argv);                               \
    std::fflush(stdout);                                                       \
    std::_Exit(rc);                                                            \
  }
#else
#define KF_TEST_MAIN()                                                         \
  int main(int argc, char **argv) { return kftest::runTests(argc, argv); }
#endif

#define KF_CHECK(cond)                                                         \
  do {                                                                         \
    if (!(cond))                                                               \
      kftest::reportFailure(__FILE__, __LINE__, "CHECK(" #cond ")");           \
  } while (0)

#define KF_CHECK_NEAR(a, b, tol)                                               \
  do {                                                                         \
    const double kf_a = (a), kf_b = (b), kf_tol = (tol);                       \
    if (!(std::abs(kf_a - kf_b) <= kf_tol)) {                                  \
      std::ostringstream kf_os;                                                \
      kf_os << #a " = " << kf_a << ", " #b " = " << kf_b << ", |diff| = "      \
            << std::abs(kf_a - kf_b) << " > " << kf_tol;                       \
      kftest::reportFailure(__FILE__, __LINE__, kf_os.str());                  \
    }                                                                          \
  } while (0)

/** Relative comparison: |a-b| <= rel * max(|a|,|b|, floor). */
#define KF_CHECK_REL(a, b, rel, floor)                                         \
  KF_CHECK_NEAR(a, b,                                                          \
                (rel) * std::max({std::abs(double(a)), std::abs(double(b)),    \
                                  double(floor)}))
