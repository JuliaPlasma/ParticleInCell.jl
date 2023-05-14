import ParticleInCell2
using PkgBenchmark

results = judge(ParticleInCell2, "main")

export_markdown(stdout, results)
