# Blitz to mdspan Migration Notes

This note records practical lessons from migrating a real scientific codebase
from Blitz++ arrays to the `kokkos/mdspan` reference implementation.

The short version is:

- `mdspan` is a good target for modern C++ array views.
- It is not a drop-in replacement for Blitz++.
- The safest migration is to separate ownership from indexing first, then
  replace Blitz-specific conveniences one by one.

## Why Migrate

Common reasons to move away from Blitz++:

- Blitz++ is old and often awkward to package on modern systems.
- The rest of the code may already be moving toward standard C++17 or later.
- `mdspan` gives a standard-style multidimensional view abstraction.
- `mdspan` plays more naturally with `std::vector`, HDF5, MPI buffers, and
  modern toolchains.

## What mdspan Does and Does Not Replace

`mdspan` replaces multidimensional indexing over contiguous storage.

It does not replace all of Blitz++:

- It does not own memory.
- It does not provide Blitz slicing syntax.
- It does not provide Blitz expression-template assignment style.
- It does not provide helper types like `TinyVector`.

That means a direct search-and-replace from `blitz::Array` to `mdspan` is not
realistic for most codes.

## Recommended Migration Strategy

Use a staged migration.

### 1. Add mdspan without changing numerics

- Vendor `kokkos/mdspan` into the repository.
- Raise the C++ baseline if needed; C++17 is a reasonable floor for the
  reference implementation.
- Keep the physics and algorithms unchanged at first.

### 2. Introduce a project-owned array wrapper

Before removing Blitz++, create your own small array type or wrapper layer.

In this repo, the critical move was replacing direct `blitz::Array` aliases
with a project wrapper in `include/decls.hpp`, backed by:

- `std::vector<T>` for ownership
- `std::mdspan` for indexing

This wrapper handled the minimum operations the code already expected:

- `resize(...)`
- `operator()(...)`
- `data()`
- `extent(int)`
- `size()`
- copy/move semantics
- legacy compatibility helpers like `lbound()` and `ubound()`

This is what made the migration tractable.

### 3. Replace helper types separately

Replace Blitz helper types independently of the array migration.

Typical replacements:

- `blitz::TinyVector<T, N>` -> `std::array<T, N>`
- Blitz dimension/bounds utilities -> project helpers or simple arrays

This was easy in comparison with the array replacement.

### 4. Rewrite slice-heavy code last

The hardest code is usually not the core storage itself. It is the Blitz slice
syntax around boundaries, ghost cells, and stencil updates.

Example pattern that does not translate directly:

```cpp
x(II, II, R1, Ic) = x(II, II, R2, Ic);
```

There is no direct `mdspan` equivalent for that style. The realistic
replacement is explicit loops or small helper routines.

Do this last, after the wrapper and indexing model are stable.

## What Was Hard in Practice

These were the hardest parts of the real migration.

### 1. Boundary-condition code using Blitz slicing

This was the biggest source of friction.

Why it was hard:

- Blitz makes ghost-cell copying and reflection compact.
- `mdspan` only gives indexed access, not slice assignment.
- Boundary code often contains off-by-one sensitive physics logic.

What worked:

- Replace each slice assignment with explicit nested loops.
- Keep the original loop bounds and component selection as intact as possible.
- Avoid combining logic cleanups with the slice rewrite.

Lesson:

- Treat boundary-condition code as high-risk.
- Expect this to be the slowest and most error-prone part of the migration.

### 2. Ownership vs view semantics

Blitz arrays often behave like a combined owner/view object. `mdspan` does not.

Why it was hard:

- Existing code expects resizeable owning arrays.
- Existing code also expects pointer-style access for HDF5 and solver kernels.
- Copying and moving arrays must leave the internal view valid.

What worked:

- Store extents separately.
- Own storage in `std::vector<T>`.
- Rebuild the `mdspan` view any time storage or shape changes.
- Implement explicit copy and move constructors/assignment operators.

Lesson:

- If you skip the ownership layer, you will fight lifetime bugs everywhere.

### 3. Preserving legacy call sites

A full API redesign during migration sounds attractive, but it multiplies risk.

Why it was hard:

- Many call sites already relied on old array behaviors.
- A cleanup-minded rewrite can silently change indexing assumptions.

What worked:

- Keep the wrapper API close to the old one first.
- Add only the minimum compatibility helpers needed to keep the code stable.
- Defer API cleanup until after regression tests pass.

Lesson:

- Compatibility glue is worth it during the migration.
- Cleanup should be a follow-up task, not part of the core port.

### 4. Testing assumptions were wrong

One migration surprise was that a seemingly convenient small regression case was
not actually a good physics baseline.

What happened here:

- A tiny `riemann1d` setup looked like a good candidate for snapshot testing.
- On closer inspection, that small case already produced NaNs even before the
  mdspan migration.
- That made it a bad regression baseline for validating the port.

What worked:

- Replace the test with a small `diffusion` case that stayed finite and
  deterministic in both serial and MPI runs.

Lesson:

- Do not assume a "small" case is a valid regression case.
- Validate the test problem itself before trusting it as a migration guard.

## Suggested Replacement Patterns

### Owning array + mdspan view

General shape:

```cpp
template <class T, std::size_t Rank>
class Array
{
public:
  void resize(...);
  T& operator()(...);
  const T& operator()(... ) const;
  T* data();
  int extent(int r) const;

private:
  std::array<int, Rank> extents_;
  std::vector<T> storage_;
  std::mdspan<T, std::dextents<int, Rank>, std::layout_right> view_;
};
```

This is the core pattern that made the migration work.

### TinyVector replacement

```cpp
std::array<int, 3> lower;
std::array<int, 3> upper;
```

### Slice replacement

Blitz style:

```cpp
x(II, II, R1, Ic) = x(II, II, R2, Ic);
```

Replacement style:

```cpp
for(int iz=0; iz < nz ;iz++) {
  for(int iy=0; iy < ny ;iy++) {
    for(int ix=0; ix < nb ;ix++) {
      x(iz, iy, dst0 + ix, ic) = x(iz, iy, src0 + ix, ic);
    }
  }
}
```

This is verbose, but it is explicit and debuggable.

## Testing Strategy

Use tests as migration gates, not as an afterthought.

Recommended order:

1. Add one stable single-process regression snapshot.
2. Add one MPI regression for the same physics case.
3. Keep one smoke test for a second executable.
4. Run tests after every array, boundary, MPI, or HDF5 change.

Good regression checks for this class of code:

- final HDF5 datasets compared against a known snapshot
- serial and MPI decomposition consistency
- stable finite outputs for representative small cases

## Practical Advice for Other Codes

If another code still uses Blitz++, a good migration checklist is:

- inventory where Blitz is used:
  - ownership
  - indexing
  - slices
  - helper types
- vendor `mdspan` first
- add a wrapper layer before replacing arrays directly
- replace `TinyVector` and similar helpers early
- rewrite slice-heavy kernels cautiously and separately
- add regression tests before changing numerics-sensitive code

## What I Would Do Earlier Next Time

- I would validate regression cases earlier, before trusting them.
- I would isolate boundary-condition rewrites into even smaller commits.
- I would keep the first wrapper even more conservative and postpone cleanup.
- I would document ownership/view assumptions at the same time as the wrapper.

## Bottom Line

The migration is very doable, but the difficulty is not in `mdspan` itself.
The difficulty is in replacing the extra conveniences Blitz++ provided around
it: ownership, helper types, and especially slicing.

If you separate those concerns and keep good regression tests in place, the
migration is manageable and reusable across similar scientific codes.
