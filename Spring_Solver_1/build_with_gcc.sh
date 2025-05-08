#!/usr/bin/bash
export PATH="/mingw64/bin:$PATH"

echo ">>> Current dir: $(pwd)"
echo ">>> Running build_with_gcc.sh"
echo ">>> Raw args: $@"

SRC_WIN=""
OUT_WIN=""

# Parse args to extract source file and output file
i=0
while [[ $i -lt $# ]]; do
  arg="${@:$((i+1)):1}"
  next="${@:$((i+2)):1}"

  echo "arg[$i] = $arg"

  if [[ "$arg" == *.c ]]; then
    SRC_WIN="$arg"
  elif [[ "$arg" == "-o" && -n "$next" ]]; then
    OUT_WIN="$next"
    i=$((i + 1))
  fi
  i=$((i + 1))
done

if [[ -z "$SRC_WIN" ]]; then
  echo "ERROR: No source file (.c) found!"
  exit 1
fi

# Convert input and output to Windows-native format for gcc
SRC_WIN_NATIVE=$(cygpath -w "$SRC_WIN") || exit 1
echo ">>> Converted source path: $SRC_WIN_NATIVE"

if [[ -n "$OUT_WIN" ]]; then
  OUT_WIN_NATIVE=$(cygpath -w "$OUT_WIN") || exit 1
  echo ">>> Converted output path: $OUT_WIN_NATIVE"
fi

# Rebuild argument list with corrected paths
args=()
for a in "$@"; do
  if [[ "$a" == "$SRC_WIN" ]]; then
    args+=("$SRC_WIN_NATIVE")
  elif [[ "$a" == "$OUT_WIN" ]]; then
    args+=("$OUT_WIN_NATIVE")
  else
    args+=("$a")
  fi
done

# Add static LAPACK, BLAS, and GFortran runtime libraries
args+=("/mingw64/lib/liblapack.a")
args+=("/mingw64/lib/libblas.a")
args+=("-lgfortran")
args+=("-lm")

echo ">>> Final gcc command:"
printf 'gcc '
for a in "${args[@]}"; do printf '"%s" ' "$a"; done
echo

# Run the compiler
gcc "${args[@]}"
