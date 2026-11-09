# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Header-only C++ library of speed-optimized DSP building blocks using x86 SIMD intrinsics. Everything lives in `include/`; there is no build system, library target, or CI. Active development is almost entirely in `include/evm.h` (efficient vector math), `include/evm_trig.h` (the SIMD sine, cosine and atan2 its polar/cartesian conversions use) and their tests. Other headers will be updated
or added over time but may not work right now.

## Building and testing

Tests are compiled by hand from inside `unit_test/` (each `.cc` has its compile line in its header comment). C++17 is required by `evm.h`.

```sh
cd unit_test
g++ -std=c++17 -march=native -O3 -I../include/ -o evm_unit_test evm_unit_test.cc
./evm_unit_test            # random seed, printed at start
./evm_unit_test 12345      # reproduce a run with a specific seed
```

Exit status is 0 only if every check passes. There is no per-test filter: the binary always runs the full suite (about 100 checks, a few seconds).

The SIMD branches are independent implementations, so test each one by forcing dispatch down with macros:

```sh
g++ ... -DDISABLE_AVX512 ...   # AVX2/AVX branch
g++ ... -DDISABLE_AVX2 ...     # AVX branch
g++ ... -DDISABLE_SIMD ...     # scalar "default" branch
g++ -std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -I../include/ -o ut_san evm_unit_test.cc
```

The untracked `ut_AVX`, `ut_AVX2`, `ut_scalar`, `ut_san` binaries in `unit_test/` are these variants. The code also targets Intel `icpx` (`-xHost` in place of `-march=native`). `evm_tail_speed_comparison.cc` is a benchmark, not a test. There is no `.gitignore`, so never stage built binaries.

Build and test with both g++ and icpx, make sure to test AVX512, AVX2, AVX, and
default paths through use of DISABLE_<ISA>. With icpx, also build with `-DDISABLE_SVML`, which switches the polar/cartesian functions from SVML to `evm_trig.h`.

## evm.h architecture

**Runtime dispatch via function multiversioning.** Each public function is defined several times with the same signature, each with `__attribute__((__target__("avx512f")))`, `"avx"` (sometimes `"avx2"` / `"avx2,fma"`), and `"default"`. GCC/icpx pick the best variant at runtime from the CPU, regardless of `-march`. Every non-default variant is wrapped in `#if !defined(DISABLE_<ISA>)`. SSE variants exist only for `getOptLevelRaw()`. Adding a function means writing the AVX512, AVX, and default variants, each with its guard.

A plain (non-multiversioned) wrapper that calls a multiversioned function, such as the by-value `dotp`/`dotpc` wrappers, `magSq(float)` and `abs`, must be defined after every target version of that function. Declare it in the top half, define it after the `default` version. Defined earlier, icpx refuses to compile the header, and g++ binds the call to the first clone (`.avx512f`), so there is no runtime dispatch.

**File layout.** The top half of the file has doc-commented declarations grouped by `//===// Section //===//` banners. The bottom half has matching `// <Section> Function Implementation //` banners with the definitions. Keep the two in sync. The file-header comment lists the operations and their accuracy and speed tradeoffs. The specialized functions (tune, fs2Shift) are listed there but not implemented yet.

**Tail handling is where bugs live.** Main loops process full registers. The final partial register uses masked loads and stores: `masks[rem]` with `_mm256_maskload_ps`/`maskstore` for AVX, and `mask16(rem)` with `maskz_loadu`/`mask_storeu` for AVX512. Complex data uses `masks[rem<<1]`. Every function must:
- treat `len <= 0` as a no-op,
- never write past `out[len-1]`,
- work in place (`out` may alias an input).

`evm_unit_test.cc` sweeps every remainder for 4/8/16-element registers, plus negative lengths, to catch these bugs. Register new functions there via `check<>()` with the right tolerance class: `TOL_EXACT` for true multiply/divide, `TOL_NR` for reciprocal plus one Newton-Raphson step, `TOL_RCP` for a bare reciprocal, `TOL_CMP` (exact) for comparison and selection. The doubles baseline is computed with `std::complex<double>`. The dot-products are reductions and go through `checkDot()` instead; they must set the result to 0 when `len <= 0`. Output buffers are pre-filled with a NaN sentinel that is matched by its bit pattern, so a NaN computed by a kernel counts as a wrong value or a stray write, not as an untouched element. The complex magnitudes run on `gZc`, which has a zero in every third element, because `x*rsqrt(x)` is `0*inf = NaN` at zero unless the kernel guards it. Its zeros cycle through all four sign combinations, because `atan2(±0, -0)` is ±π (`cartToPol` runs on it too). The comparison functions run on `gNf`, whose values are a multiple of 0.5 half of the time, so ties with the threshold and equal magnitudes of opposite sign are common. `minMax` with two output buffers goes through `checkPair()`, which also runs each output on top of either input, and the `minMax` reduction goes through `checkMinMax()`, which puts a single min or max at every position. Where a stray write could put back the sentinel unchanged (min/max of two sentinels, `zneg`/`zpos` in place), those tests fill the area past the end with values the kernel would change. `multSum` and `multcSum` also pass `check<>()` a magnitude function, |in1*in2|+|in3|, to measure the error against, because their results can cancel to near zero; their "a=in3" runs pass in3 as the a buffer so the in-place modes also cover out on top of in3. `multe`, `dive`, `divre` and `divnre` go through `checkExpand()`, which requires every imaginary part to be left exactly as it was (the sentinel, or in place the input float under it) and runs with the complex output over in1 or in2. The divide-gather functions (`divg`, `divrg`, `divnrg`) go through `checkGather()`, which runs strides from 1 to INT_MAX on the 16K-element `gGa`/`gGb`, each output count with the shortest and the longest input that gives it, and requires everything past `out[n-1]` to be left as it was (in place, the rest of the input); `len <= 0` or `stride <= 0` must be a no-op. On unix-like systems it also runs each case with the inputs ending right before a `PROT_NONE` page and catches the SIGSEGV. The small strides load whole registers, so a read past the end would otherwise go unseen. Their eps overloads and `divSqrt`'s use `GEPS` = 0.5, which raises every negative divisor, and `divSqrt` without eps runs on `gPf` = |`gBf`|. The integer to float conversions run on `gI8`/`gI16`/`gI32` and must equal `static_cast<float>` exactly; their "in place" runs go through the `void*` wrappers, whose float output is wider than the integer input, and the min/max overloads go through `checkCnvrtMinMax()`. `polToCart` always uses `evm_trig.h`; the other polar/cartesian SIMD versions use SVML when `EVM_SVML` is defined (icpx without `-DDISABLE_SVML`) and `evm_trig.h` otherwise. Their error is measured against `TOL_TRIG`. The polar ones also run on large phases (`gLc` up to 1e4 cycles, `gLr` up to 1e7 radians, far past the 65536 where `evm_trig.h` hands lanes over to scalar `sincosf`), and `checkPolarSpecial()` runs grids of zeros of both signs, ±1, infinities, NaN, denormals and the largest floats through all four, requiring NaN and infinities exactly and the sign of a zero angle. icpx turns on FTZ/DAZ by default, so there denormal inputs act as zeros (the double references see them that way too) and denormal results may come out as 0; `check<>()` accepts a 0 for a result below FLT_MIN only while flush to zero is on (it reads MXCSR), so g++ builds still check denormal results.

**Deliberate choices that look like bugs:**
- Division kernels divide the zero-filled unused lanes of the last register (0/0), which raises `FE_INVALID`/`FE_DIVBYZERO` on valid input. This is intentional for speed and results are unaffected. `setUnusedToOne()` exists for callers who need clean FP flags.
- "True" complex division computes `1/divisor` and then multiplies. This is faster but less accurate (about 1.5 ulp) and range-restricted.
- `_mm256_blendv_ps` is avoided because GCC 12 can expand it into a branch per element.
- `min_ps`/`max_ps` return their second operand on a tie or NaN, so the min/max kernels pass `in2` first (and `zneg`/`zpos` pass 0 first) to match `std::min(in1, in2)` in the default versions. The `minMax` reduction passes the loaded value first so that a NaN element is skipped.
- `multSum` and `multcSum` results can differ in the last bit between dispatch paths: the AVX-512 and AVX2 versions use FMA (one rounding), the AVX and default versions round the product first. Built with `-march=native`, GCC can fuse the AVX version too.
- The int8/int16 to float conversions run from the end of the buffer to the start, so the in-place `void*` wrappers work even though the output is wider than the input. Their AVX512 versions also require AVX512BW and AVX512VL, for masked byte/word loads and 512-bit byte/word min/max (the min/max loops cover 64 or 32 elements per pass): the target is `"avx512f,avx512bw,avx512vl"`, since g++ 12 rejects `"avx512bw"` on its own in a dispatched target. AVX and AVX2 have no masked byte/word loads, so their versions do the partial register first, by loading the last full register (overlapping the one before) and storing only its last `rem` floats with `rmasks[rem]`, and convert buffers shorter than one register one element at a time. Their AVX2 min/max versions also find the min/max on the raw bytes/words (256-bit, 32 or 16 elements per pass); each remainder register is loaded for the min/max with a broadcast that repeats its 8 integers across the register, and the main pass loads and converts every register before storing any of them, which is 10-17% faster with icpx than alternating. The scalar int16 loops read the input with `memcpy`: in place, the `int16_t` input and `float` output alias, and otherwise icpx assumes they don't and reorders the loop.
- `multe`, `dive`, `divre` and `divnre` run from the end of the buffer to the start, so they work in place with the complex output over in1 or in2, and their SIMD versions store only the real parts (AVX-512 mask `0x5555`, AVX/AVX2 `maskstore` of every other float), so the imaginary parts are never written. On TGL the masked stores were also faster than loading the output and blending the products in.
- `cnvrtInt32ToFloat` without min/max, `remult` and `immult` have no AVX2 version: their AVX versions already use only the instructions an AVX2 one would. `dive` and `divSqrt` have none either: the divider limits them, so AVX2's permutes gained little over AVX's unpack and permute2f128 (on TGL, `dive` 0.62 vs 0.63 cycles per element, `divSqrt` 1.37 vs 1.39-1.41). The polar/cartesian AVX2 versions are the AVX ones recompiled (plus FMA in the `cartToPol` magnitude), because SVML's AVX2 `sincos` is about 2.4x and its `atan2` 1.2-1.4x faster than its AVX ones.
- On the SVML path, `cartToPol`/`cartToPolRad` set re to ±1 wherever `re|im` compares equal to 0 (the bits of 1.0 or'd into re) before the `atan2`, which gives the same result, and their default versions move re away from 0 by ±1 wherever im is 0 (branch-free, since the angle only depends on re's sign there), since icpx vectorizes those loops with SVML. SVML's AVX-512 `atan2` takes a slow path for any register with a (0,0) input (1.8x slower with one per 37 elements, 5.3x on zero-padded data), and every SVML `atan2` returns NaN for a both-denormal pair when denormals are treated as zero (DAZ, icpx's default). The float compare, rather than a bit test, is what also catches the denormal pairs under DAZ. Without DAZ, both-denormal inputs still hit SVML's slow path; they are left alone, since their angle is real. `evm_trig.h`'s `atan2` has neither problem.
- `evm_trig.h` uses the Cephes single precision polynomials. `sincosCyc` reduces the phase exactly (to the nearest quarter cycle), so `polToCart` keeps full accuracy at any phase, unlike SVML or `sincosf(2*pi*ph)`; the default `polToCart` removes whole cycles before `sincosf` for the same reason. Its AVX2/AVX versions give NaN past |ph| = 8.5e37, where `4*ph` overflows (removing whole cycles first would cost 13-20%); the AVX-512 one rounds with `roundscale` and has no limit. `sincosRad` reduces with pi/2 in three pieces with FMA (five 8-bit pieces on AVX) and is accurate to |x| = 65536, beyond which those lanes fall back to scalar `sincosf`. Accuracy is about 1.7 ulp for sin/cos (3.3 for radians on AVX) and 2.5 ulp for `atan2`, against 0.6-0.9 for SVML, which is why the radians and `atan2` functions keep SVML when it is available.
- `multDivnr`, `divnre` and `divnrg` apply their Newton-Raphson stage to the quotient on the AVX-512 and AVX2 paths (for `multDivnr`, `q0 + r*(in2 - in3*q0)` with `q0 = in2*r`), rather than to the reciprocal as the `divnr` kernels do. With FMA it costs the same and is more accurate: on AVX-512 it equals the true divide 99.8% of the time (`divnr`'s reciprocal form about 82%). The AVX versions have no FMA, so they refine the reciprocal like `divnr`.
- `divg`, `divrg` and `divnrg` run front to back, which works in place since output k comes from input k*stride >= k. They don't use the gather instructions. With the Gather Data Sampling microcode mitigation (Intel Skylake through Ice Lake/Tiger Lake/Rocket Lake), each gather, ymm or zmm, takes about 28 cycles, which made them 2.2x (AVX-512) and 4x (AVX2) slower than scalar loads on TGL. Instead, the shared helpers `stridedLoad512`, `stridedLoadAvx2` and `stridedLoadAvx` (defined just before them) take small strides out of whole registers with permutes: AVX-512 uses a chain of `vpermt2ps` up to `STRIDED_PERM_MAX` = 8 (indices from `stridedIdx512`, set once per call; each step only sets its own register's lanes, so a `switch` with fallthrough runs them, which is 25% faster at strides 3-4 than a loop, where g++ tests the stride after every step), AVX2 uses fixed shuffles for strides 2-4, and AVX for stride 2. Other strides load each float and insert it (`_mm256_setr_ps`). At strides 2-4 this runs at the divider's speed on AVX-512 (0.64 cycles per output, 7x faster than the gathers). The main loop only takes registers whose whole input span (16 or 8 times the stride) lies inside the buffers. The last 1-16 (or 1-8) outputs use masked loads, or in the insert case lanes that repeat the last float, so nothing past the end is read. With buffers that aren't 64-byte aligned and come from L2, the permutes at strides 6-8 are up to 30% slower than inserts, since every register load then spans two cache lines (aligned, they're 25-45% faster). A stride <= 0 is a no-op in every version, including the defaults (where a stride of 0 would never end).
- `divSqrt` computes the square root and then its reciprocal with one divide, and multiplies by that, like the complex/real `div` (see "true" complex division above). `divnrSqrt` refines `rsqrt` with `r + 0.5*r*(1 - x*r*r)` on the FMA paths, which is more accurate than `0.5*r*(3 - x*r*r)` (on AVX-512, 1 ulp of 1/sqrt vs 2, and correctly rounded 87% of the time vs 75%). The AVX version uses the latter. The eps overloads, and `divSqrt`, `divrSqrt` and `divnrSqrt`, clamp the divisor with `max_ps(eps, x)` (x second, see min/max above) to match `std::max(x, eps)`.
- AVX/AVX2 `magr` computes `|in|` as `rcp(rsqrt(|in|^2))`, not `|in|^2*rsqrt(|in|^2)`. It costs the same, but it gives 0 at zero instead of `0*inf = NaN`. AVX512 uses a zero-masked `rsqrt14` instead, because `rcp14` is 3 uops on a zmm register.
- getOptLevel SSE variants are intentional even though there are no
  SSE-optimized function

**icpx compatibility constraints** (see "Compiler-dependent notes" in the header):
- `default` variants must not use the `std::complex<float>(re, im)` constructor, `conj()`, or float-times-complex arithmetic, because icpx leaves them as undefined references inside target-versioned functions. Expand with `.real()`/`.imag()` or copy a constexpr complex and set its parts.
- `std::string` can't be constructed inside multiversioned functions. This is why `getOptLevel()` wraps `getOptLevelRaw()`.
- Multiversioned functions can't have default arguments. icpx doesn't apply a default argument from the declaration to the target versions, so a call that leaves it out fails to compile (g++ accepts it). Use a plain overload without the argument, declared next to the full one and defined after its `default` version (`divSqrt`, `divrSqrt` and `divnrSqrt` without `eps`).

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
- Try to keep fixes compiler-neutral, tell the user if there specific fixes for
  GCC vs. Intel compiler but favor generic fixes instead of specific ones.
- Do not commit or stage any changes to git, the user will handle that.