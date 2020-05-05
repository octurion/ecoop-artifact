# Modifications

Our modifications on the original code consist of the following:
- We eliminated the use of a nonstandard `gcc` C++ extension (empty structs),
  which was causing us compiler errors in our version of `gcc`. The
  `p0.patch` file consists of the relevant changes (these can also be seen
  through `git diff HEAD~1`.
- We added a `script.sh` file in the root directory for automating the
  generation of the results.
