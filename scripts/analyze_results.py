import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
CSV_PATH = (BASE_DIR / "../outputs/runtime_log.csv").resolve()
PLOTS_DIR = (BASE_DIR / "../plots").resolve()
PAPER_DIR = (BASE_DIR / "../paper").resolve()

PLOTS_DIR.mkdir(parents=True, exist_ok=True)
PAPER_DIR.mkdir(parents=True, exist_ok=True)

if not CSV_PATH.exists():
    raise FileNotFoundError(f"Runtime log not found at {CSV_PATH}")

df = pd.read_csv(CSV_PATH)
df["real_time_sec"] = pd.to_numeric(df["real_time_sec"], errors="coerce")

serial_runtimes = df[df["np"] == 2].groupby("dataset")["real_time_sec"].min().to_dict()

def _compute_speedup(row):
    baseline = serial_runtimes.get(row["dataset"], row["real_time_sec"])
    if pd.isna(row["real_time_sec"]) or pd.isna(baseline) or row["real_time_sec"] == 0:
        return float("nan")
    return baseline / row["real_time_sec"]


df["speedup"] = df.apply(_compute_speedup, axis=1)
df["efficiency"] = df["speedup"] / df["np"]

# Tables
table_runtime = df.pivot_table(index="np", columns="dataset", values="real_time_sec").sort_index()
table_speedup = df.pivot_table(index="np", columns="dataset", values="speedup").sort_index()
table_eff = df.pivot_table(index="np", columns="dataset", values="efficiency").sort_index()

(PAPER_DIR / "runtime_table.md").write_text(table_runtime.round(2).to_markdown())
(PAPER_DIR / "speedup_efficiency_table.md").write_text(
    "# Speedup Table\n" + table_speedup.round(2).to_markdown() + "\n\n# Efficiency Table\n" + table_eff.round(2).to_markdown()
)

# Plots
sns.set(style="whitegrid")
metric_to_filename = {
    "real_time_sec": "runtime_plot.png",
    "speedup": "speedup_plot.png",
    "efficiency": "efficiency_plot.png",
}

for metric, filename in metric_to_filename.items():
    plt.figure(figsize=(8, 5))
    for dataset in sorted(df["dataset"].unique()):
        sub = df[df["dataset"] == dataset].sort_values("np")
        plt.plot(sub["np"], sub[metric], marker="o", label=dataset)
    plt.xlabel("Number of Processes")
    plt.ylabel(metric.replace("_", " ").title())
    plt.title(f"{metric.replace('_', ' ').title()} vs. #Processes")
    plt.legend()
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / filename, dpi=150)
    plt.close()
