# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Header-only C++ library of speed-optimized DSP building blocks using x86 SIMD intrinsics. Everything lives in `include/`; there is no build system, library target, or CI. Active development is almost entirely in `include/evm.h` (efficient vector math) and its tests. Other headers will be updated
or added over time but may not work right now.

## Building and testing

Tests are compiled by hand from inside `unit_test/` (each `.cc` has its compile line in its header comment). C++17 is required by `evm.h`.

```sh
cd unit_test
g++ -std=c++17 -march=native -O3 -I../include/ -o evm_unit_test evm_unit_test.cc
./evm_unit_test            # random seed, printed at start
./evm_unit_test 12345      # reproduce a run with a specific seed
```

Exit status is 0 only if every check passes. There is no per-test filter: the binary always runs the full suite (about 38 checks, a few seconds).

The SIMD branches are independent implementations, so test each one by forcing dispatch down with macros:

```sh
g++ ... -DDISABLE_AVX512 ...   # AVX2/AVX branch
g++ ... -DDISABLE_AVX2 ...     # AVX branch
g++ ... -DDISABLE_SIMD ...     # scalar "default" branch
g++ -std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -I../include/ -o ut_san evm_unit_test.cc
```

The untracked `ut_AVX`, `ut_AVX2`, `ut_scalar`, `ut_san` binaries in `unit_test/` are these variants. The code also targets Intel `icpx` (`-xHost` in place of `-march=native`). `evm_tail_speed_comparison.cc` is a benchmark, not a test. There is no `.gitignore`, so never stage built binaries.

Build and test with both g++ and icpx, make sure to test AVX512, AVX2, AVX, and
default paths through use of DISABLE_<ISA>.

## evm.h architecture

**Runtime dispatch via function multiversioning.** Each public function is defined several times with the same signature, each with `__attribute__((__target__("avx512f")))`, `"avx"` (sometimes `"avx2"` / `"avx2,fma"`), and `"default"`. GCC/icpx pick the best variant at runtime from the CPU, regardless of `-march`. Every non-default variant is wrapped in `#if !defined(DISABLE_<ISA>)`. SSE variants exist only for `getOptLevelRaw()`. Adding a function means writing the AVX512, AVX, and default variants, each with its guard.

**File layout.** The top half of the file has doc-commented declarations grouped by `//===// Section //===//` banners. The bottom half has matching `// <Section> Function Implementation //` banners with the definitions. Keep the two in sync. The file-header comment lists the operations and their accuracy and speed tradeoffs. Dot-product, power, magnitude, and specialized (tune, fs2Shift) are listed there but not implemented yet: their implementation sections are empty, and the Dot-Product declarations are in progress.

**Tail handling is where bugs live.** Main loops process full registers. The final partial register uses masked loads and stores: `masks[rem]` with `_mm256_maskload_ps`/`maskstore` for AVX, and `mask16(rem)` with `maskz_loadu`/`mask_storeu` for AVX512. Complex data uses `masks[rem<<1]`. Every function must:
- treat `len <= 0` as a no-op,
- never write past `out[len-1]`,
- work in place (`out` may alias an input).

`evm_unit_test.cc` sweeps every remainder for 4/8/16-element registers, plus negative lengths, to catch these bugs. Register new functions there via `check<>()` with the right tolerance class: `TOL_EXACT` for true multiply/divide, `TOL_NR` for reciprocal plus one Newton-Raphson step, `TOL_RCP` for a bare reciprocal. The doubles baseline is computed with `std::complex<double>`.

**Deliberate choices that look like bugs:**
- Division kernels divide the zero-filled unused lanes of the last register (0/0), which raises `FE_INVALID`/`FE_DIVBYZERO` on valid input. This is intentional for speed and results are unaffected. `setUnusedToOne()` exists for callers who need clean FP flags.
- "True" complex division computes `1/divisor` and then multiplies. This is faster but less accurate (about 1.5 ulp) and range-restricted.
- `_mm256_blendv_ps` is avoided because GCC 12 can expand it into a branch per element.
- getOptLevel SSE variants are intentional even though there are no
  SSE-optimized function

**icpx compatibility constraints** (see "Compiler-dependent notes" in the header):
- `default` variants must not use the `std::complex<float>(re, im)` constructor, `conj()`, or float-times-complex arithmetic, because icpx leaves them as undefined references inside target-versioned functions. Expand with `.real()`/`.imag()` or copy a constexpr complex and set its parts.
- `std::string` can't be constructed inside multiversioned functions. This is why `getOptLevel()` wraps `getOptLevelRaw()`.

**Namespace hygiene.** Keep `using std::complex` and other `using` declarations inside `namespace EVM`, never at global scope (commit f6c7ce4), so they don't leak into users of `evm.h`.

**Memory.** `aligned_allocator.h` gives STL-compatible aligned allocation. Alignment and pad are specified in bytes. `ALIGNED_VECTOR(T)` uses `evm_align`, which the `DISABLE_*` macros select: 64-byte alignment plus 2 registers of trailing pad when AVX512 is enabled. `EVM::getSize<T>(n)` rounds a length up to a whole number of registers. `EVM::enableFastDenormals()` should be called at program start.

**Processing scenarios** evm.h normally runs on longer buffers and on AVX2 or AVX-512 CPUs. When an optimization trades a few percent at short lengths (under ~32 elements) or on the AVX-only path for a long-buffer gain, the user prefers the simpler code that wins on long buffers.  Also prioritize absolute speed over FE_INVALID/FE_DIVBYZERO-clean tails.  Weight long-buffer AVX2/AVX-512 results most when recommending kernel changes; mention short-length or AVX-only regressions and let the user decide if they are worth fixing.  The exception is
dot-product methods which will normally be run on shorter buffers.

**Performance benchmark pitfalls**
- 4K aliasing: separate large std::vectors land at the same offset within a 4 KB page, so back-to-back calls serialize (next call's loads wait on the previous call's stores) and longer-latency kernels look slower at 8–64 elements. Time both that layout and buffers offset ~1 KB apart in one pool.
- Code placement: differences under ~0.3 ns per call (about 1 cycle) can come from where code lands; identical machine code measured 3.24 vs 3.46 ns in two builds. Compare versions in the same binary, best of 3+, and don't attribute such gaps to source changes without checking the disassembly.
- Extracted test copies need their `default` version too, or GCC compiles them as plain functions called directly (no ifunc dispatch), which flatters short-call timings


## Other headers

These predate evm.h and are not exercised by the current tests:
- `genfilter.h` (namespace `GenFilter`): windowed-sinc FIR design.
- `polyphase.h`: polyphase filter bank from such filters.
- `resampler.h` / `resamplerps.h`: arbitrary-rate resampling. Both define `class Resampler`, so never include both in one translation unit. They still refer to old `-DHAVE_SSE*` / `HAVE_CONFIG_H` macros and an `inc/config.h` that doesn't exist.
- `constel.h` (namespace `Constel`): PSK/QAM constellations and bit mapping.
- `timecode.h`: sample-indexed timestamps.
- `equalizer.h` (LMS) and `dft_wrapper.h`: depend on Intel MKL (`mkl.h`).

`dft_wrapper.h` is a non-compiling stub, and `unit_test/constel_test.cc` is unfinished and does not compile.

## Reviewing Diffs

Present proposed code changes as a graphical VS Code diff: `code --reuse-window --diff <read-only snapshot of the current file> <my modified copy in the scratchpad>`. Left = the old version, right = my modified version (editable). The user rejects any change by clicking that block's arrow (VS Code's "Revert Block" copies the left text back in to the right file), keeps the rest, saves the right file and tells me. I then check the real file is unchanged since the snapshot (save its checksum when making the snapshot), copy the reviewed file into the project, and rebuild/retest.

**Why:** the user wants to see old vs new and revert individual changes.

**How to apply:** prepare and test the change in a scratch copy first, chmod the snapshot read-only, and when opening the diff say plainly that the arrow rejects a change. Verify the reviewed file against my proposal before installing it and report any blocks the user reverted.

## Conventions

- Each source file starts with the LGPL v3 license block (`Copyright (c) <year> Nick Xenias`), followed by a `//===` banner with Name / Purpose / Created / Description.
- Commit messages are prefixed `Mod:`, `New:`, `Fix:`, or `Doc:`.
- Style: 2-space indent, lines kept to about 80 columns, `T const * const` for input pointers and `T * const` for outputs.
- keep fixes compiler-neutral, do not use GCC-only asm pinning
