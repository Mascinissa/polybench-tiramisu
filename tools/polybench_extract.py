"""Shared parsing utilities for PolyBench/C 4.2.1 sources and the
polybench_tiramisu generator files.

Used by:
- generate_polybench_init_sidecars.py (sidecar generation)
- validate_init_parity.py             (numeric parity vs. stock PolyBench)
- validate_against_polybench_c.py     (end-to-end kernel output check)
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

# tiramisu benchmark name -> path of the stock benchmark (relative to the
# PolyBenchC-4.2.1 root, without extension; <path>.c and <path>.h exist).
BENCHMARKS = {
    "2mm": "linear-algebra/kernels/2mm/2mm",
    "3mm": "linear-algebra/kernels/3mm/3mm",
    "adi": "stencils/adi/adi",
    "atax": "linear-algebra/kernels/atax/atax",
    "bicg": "linear-algebra/kernels/bicg/bicg",
    "cholesky": "linear-algebra/solvers/cholesky/cholesky",
    "correlation": "datamining/correlation/correlation",
    "covariance": "datamining/covariance/covariance",
    "deriche": "medley/deriche/deriche",
    "doitgen": "linear-algebra/kernels/doitgen/doitgen",
    "durbin": "linear-algebra/solvers/durbin/durbin",
    "fdtd_2d": "stencils/fdtd-2d/fdtd-2d",
    "floyd_warshall": "medley/floyd-warshall/floyd-warshall",
    "gemm": "linear-algebra/blas/gemm/gemm",
    "gemver": "linear-algebra/blas/gemver/gemver",
    "gesummv": "linear-algebra/blas/gesummv/gesummv",
    "gramschmidt": "linear-algebra/solvers/gramschmidt/gramschmidt",
    "heat3d": "stencils/heat-3d/heat-3d",
    "jacobi1d": "stencils/jacobi-1d/jacobi-1d",
    "jacobi2d": "stencils/jacobi-2d/jacobi-2d",
    "lu": "linear-algebra/solvers/lu/lu",
    "ludcmp": "linear-algebra/solvers/ludcmp/ludcmp",
    "mvt": "linear-algebra/kernels/mvt/mvt",
    "nussinov": "medley/nussinov/nussinov",
    "seidel2d": "stencils/seidel-2d/seidel-2d",
    "symm": "linear-algebra/blas/symm/symm",
    "syr2k": "linear-algebra/blas/syr2k/syr2k",
    "syrk": "linear-algebra/blas/syrk/syrk",
    "trisolv": "linear-algebra/solvers/trisolv/trisolv",
    "trmm": "linear-algebra/blas/trmm/trmm",
}

# tiramisu dataset suffix -> PolyBench dataset macro prefix
DATASETS = {
    "MINI": "MINI",
    "SMALL": "SMALL",
    "MEDIUM": "MEDIUM",
    "LARGE": "LARGE",
    "XLARGE": "EXTRALARGE",
}

# tiramisu buffer base name (b_ prefix stripped) -> stock PolyBench array
# name, where they differ. Everything else maps by identity.
BUFFER_ALIASES: dict[str, dict[str, str]] = {
    "3mm": {"E": "G"},  # tiramisu's E (NIxNL output) is stock G; stock E/F
    # are the intermediate products (tiramisu b_AB/b_CD temporaries).
    "doitgen": {"x": "C4"},
    "fdtd_2d": {"fict": "_fict_"},
    "floyd_warshall": {"paths": "path"},
    "mvt": {"y1": "y_1", "y2": "y_2"},
}

P_TYPE_TO_CTYPE = {
    "p_float64": "double",
    "p_float32": "float",
    "p_int32": "int",
}

DATA_TYPE_INFO = {
    "DOUBLE": {
        "ctype": "double",
        "printf": '"%0.2lf "',
        "scalar_val": "x",
        "sqrt": "sqrt(x)",
        "exp": "exp(x)",
        "pow": "pow(x,y)",
    },
    "FLOAT": {
        "ctype": "float",
        "printf": '"%0.2f "',
        "scalar_val": "x##f",
        "sqrt": "sqrtf(x)",
        "exp": "expf(x)",
        "pow": "powf(x,y)",
    },
    "INT": {
        "ctype": "int",
        "printf": '"%d "',
        "scalar_val": "x",
        "sqrt": "sqrt(x)",
        "exp": "exp(x)",
        "pow": "pow(x,y)",
    },
}


@dataclass
class ArrayParam:
    """An array parameter of init_array/print_array."""

    elem_type: str  # type token as written: DATA_TYPE or a typedef (base)
    name: str  # stock array name
    macro_dims: list[str]  # uppercase dim macros, e.g. ["NI", "NJ"]
    lower_dims: list[str]  # lowercase runtime dims, e.g. ["ni", "nj"]


@dataclass
class FuncInfo:
    """Parsed init_array or print_array."""

    int_params: list[str] = field(default_factory=list)  # e.g. ["ni", "nj"]
    scalar_params: list[str] = field(default_factory=list)  # e.g. ["alpha"]
    arrays: list[ArrayParam] = field(default_factory=list)
    body: str = ""  # verbatim body, including outer braces


@dataclass
class StockBenchmark:
    name: str
    c_path: Path
    h_path: Path
    data_type: str  # DOUBLE / FLOAT / INT
    typedefs: dict[str, str]  # e.g. {"base": "char"}
    datasets: dict[str, dict[str, int]]  # dataset -> {MACRO: value}
    init: FuncInfo = None  # type: ignore[assignment]
    print_: FuncInfo = None  # type: ignore[assignment]


@dataclass
class TiramisuVariant:
    benchmark: str  # e.g. "gemm"
    dataset: str  # tiramisu suffix, e.g. "MINI"
    function_name: str  # e.g. "function_gemm_MINI"
    generator_path: Path
    buffer_names: list[str]  # codegen order, e.g. ["b_A", "b_B", "b_C"]
    buffer_dims: list[list[int]]
    buffer_ctypes: list[str]


def _matched_span(text: str, start: int, open_ch: str, close_ch: str) -> tuple[int, int]:
    """Return (start, end) span of the balanced group starting at
    text[start] == open_ch; end is the index just past the closing char."""
    assert text[start] == open_ch, (text[start], open_ch)
    depth = 0
    for i in range(start, len(text)):
        c = text[i]
        if c == open_ch:
            depth += 1
        elif c == close_ch:
            depth -= 1
            if depth == 0:
                return start, i + 1
    raise ValueError(f"Unbalanced {open_ch}...{close_ch} starting at {start}")


def _split_top_level_commas(s: str) -> list[str]:
    parts, depth, cur = [], 0, []
    for c in s:
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
        if c == "," and depth == 0:
            parts.append("".join(cur))
            cur = []
        else:
            cur.append(c)
    if cur:
        parts.append("".join(cur))
    return [p.strip() for p in parts if p.strip()]


def _parse_function(source: str, func_name: str) -> FuncInfo:
    m = re.search(r"void\s+" + func_name + r"\s*\(", source)
    if not m:
        raise ValueError(f"{func_name} not found")
    paren_start = source.index("(", m.start())
    p0, p1 = _matched_span(source, paren_start, "(", ")")
    params_str = source[p0 + 1 : p1 - 1]
    brace_start = source.index("{", p1)
    b0, b1 = _matched_span(source, brace_start, "{", "}")
    body = source[b0:b1]

    info = FuncInfo(body=body)
    for param in _split_top_level_commas(params_str):
        param = " ".join(param.split())  # normalize whitespace
        m_int = re.fullmatch(r"int\s+(\w+)", param)
        if m_int:
            info.int_params.append(m_int.group(1))
            continue
        m_scalar = re.fullmatch(r"(\w+)\s*\*\s*(\w+)", param)
        if m_scalar:
            info.scalar_params.append(m_scalar.group(2))
            continue
        m_arr = re.fullmatch(
            r"(\w+)\s+POLYBENCH_([1-5])D\s*\(([^)]*)\)", param
        )
        if m_arr:
            elem_type = m_arr.group(1)
            k = int(m_arr.group(2))
            args = [a.strip() for a in m_arr.group(3).split(",")]
            # POLYBENCH_kD(name, D1..Dk, d1..dk)
            if len(args) != 1 + 2 * k:
                raise ValueError(f"Unexpected POLYBENCH_{k}D args: {args}")
            info.arrays.append(
                ArrayParam(
                    elem_type=elem_type,
                    name=args[0],
                    macro_dims=args[1 : 1 + k],
                    lower_dims=args[1 + k :],
                )
            )
            continue
        raise ValueError(f"Unrecognized {func_name} parameter: {param!r}")
    return info


def parse_stock_benchmark(polybench_root: Path, name: str) -> StockBenchmark:
    rel = BENCHMARKS[name]
    c_path = polybench_root / (rel + ".c")
    h_path = polybench_root / (rel + ".h")
    c_src = c_path.read_text()
    h_src = h_path.read_text()

    # Default data type: the header defines DATA_TYPE_IS_<T> when none of
    # the DATA_TYPE_IS_* overrides is given on the command line.
    m = re.search(r"#\s*define\s+DATA_TYPE_IS_(\w+)", h_src)
    if not m:
        raise ValueError(f"{h_path}: DATA_TYPE_IS_* not found")
    data_type = m.group(1)
    if data_type not in DATA_TYPE_INFO:
        raise ValueError(f"{h_path}: unsupported DATA_TYPE_IS_{data_type}")

    # Local typedefs used in signatures (e.g. nussinov's `typedef char base;`)
    typedefs = dict(re.findall(r"typedef\s+(\w+)\s+(\w+)\s*;", c_src))
    typedefs = {alias: underlying for underlying, alias in typedefs.items()}

    # Dataset size macros.
    datasets: dict[str, dict[str, int]] = {}
    for block_m in re.finditer(
        r"#\s*ifdef\s+(\w+)_DATASET(.*?)#\s*endif", h_src, re.DOTALL
    ):
        ds = block_m.group(1)
        macros = {
            mm.group(1): int(mm.group(2))
            for mm in re.finditer(
                r"#\s*define\s+(\w+)\s+(\d+)", block_m.group(2)
            )
        }
        if macros:
            datasets[ds] = macros

    bench = StockBenchmark(
        name=name,
        c_path=c_path,
        h_path=h_path,
        data_type=data_type,
        typedefs=typedefs,
        datasets=datasets,
    )
    bench.init = _parse_function(c_src, "init_array")
    bench.print_ = _parse_function(c_src, "print_array")
    return bench


def parse_tiramisu_variant(generator_path: Path) -> TiramisuVariant:
    src = generator_path.read_text()
    function_name = re.findall(r'tiramisu::init\("(\w+)"\)', src)[0]
    m = re.fullmatch(r"function_(\w+)_([A-Z]+)", function_name)
    if not m or m.group(2) not in DATASETS:
        raise ValueError(f"Unexpected function name: {function_name}")
    benchmark, dataset = m.group(1), m.group(2)

    codegen_line = re.findall(r"tiramisu::codegen\(\{(.+?)\}", src)[0]
    buffer_names = re.findall(r"\w+", codegen_line)

    buffer_dims: list[list[int]] = []
    buffer_ctypes: list[str] = []
    for buf in buffer_names:
        decl = re.findall(r"buffer\s+" + buf + r"\s*\(.*", src)[0]
        sizes = re.findall(r"\{(.*?)\}", decl)[0]
        buffer_dims.append([int(x) for x in re.findall(r"\d+", sizes)])
        p_type = re.findall(r"(p_\w+)", decl)[0]
        if p_type not in P_TYPE_TO_CTYPE:
            raise ValueError(f"{generator_path}: unsupported type {p_type}")
        buffer_ctypes.append(P_TYPE_TO_CTYPE[p_type])

    return TiramisuVariant(
        benchmark=benchmark,
        dataset=dataset,
        function_name=function_name,
        generator_path=generator_path,
        buffer_names=buffer_names,
        buffer_dims=buffer_dims,
        buffer_ctypes=buffer_ctypes,
    )


def stock_array_name(benchmark: str, buffer_name: str) -> str:
    """Map a tiramisu buffer name (b_X) to the stock PolyBench array name."""
    base = buffer_name[2:] if buffer_name.startswith("b_") else buffer_name
    return BUFFER_ALIASES.get(benchmark, {}).get(base, base)


def array_written_in_body(body: str, array_name: str) -> bool:
    """True if `array_name` is assigned to in the function body (heuristic:
    the identifier appears at all — PolyBench init/print bodies only
    mention arrays they access)."""
    return re.search(r"\b" + re.escape(array_name) + r"\b", body) is not None
