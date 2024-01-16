import ParticleInCell
using PkgBenchmark

results = judge(ParticleInCell, "main")

export_markdown(stdout, results)
