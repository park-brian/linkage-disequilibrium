# Profiling scripts (read-only)

Used to build [`docs/DATA.md`](../../docs/DATA.md) and [`docs/REQUIREMENTS.md`](../../docs/REQUIREMENTS.md).
Run them from `D:/Projects/LDlink` (paths are relative to the workspace root). Outputs are in `results/`.

| script | what it does |
|---|---|
| `mongo_probe.py` | Inside an ldlink api container: counts, indexes, samples, `_id` time ranges per collection. Never prints credentials. (Live hosts no longer resolve; kept for reference.) |
| `bson_profile.py` | Streams a mongodump `.bson.gz` from stdin: doc count, field types, `"NA"` counts, numeric ranges, low-cardinality values. `aws s3 cp s3://nci-cbiit-dceg-data-nonprod/ldlink/deploy/mongodb/<c>.bson.gz - \| python -I bson_profile.py <c>` |
| `vcf_profile.py` | Header/sample order vs `ALL.txt`, record classes, duplicate positions, GT forms. `python -I vcf_profile.py <vcf.gz> db/data/1000G_population_samples/ALL.txt [max_records]` |
| `geno_layout_estimate.py` | Bytes/variant for packed-bitset vs sparse-index haplotype encodings (raw, deflate, zstd). |
| `ldmatrix_reference.js` | Line-for-line port of LDlink's LDmatrix math over a local tabix-indexed VCF → reference test vectors. `node ldmatrix_reference.js convert-vcf-json/node_modules <vcf.gz> <pop.txt>... -- chr:pos:rsid ...` (`TBI=` env to point at a separate index). |
| `gt_compare_builds.js` | Per-sample GT agreement between GRCh37 phase 3 and GRCh38 High Coverage at a few positions. |

Temporary-ECS-task recipe used for the Mongo probe (always stop it afterwards):
`aws ecs run-task --cluster analysistools-dev --task-definition dev-ldlink-api --launch-type FARGATE --enable-execute-command --started-by ldlink-mongo-probe --network-configuration "awsvpcConfiguration={subnets=[subnet-409a0424,subnet-4f35f112],securityGroups=[sg-0af8bf5fca57a6be2],assignPublicIp=DISABLED}" --overrides '{"containerOverrides":[{"name":"api","command":["sleep","3600"]}]}'`
→ `aws ecs execute-command … --container api --interactive --command …` (use `python3.13`, which has pymongo)
→ `aws ecs stop-task …` and confirm `list-tasks --started-by ldlink-mongo-probe --desired-status RUNNING` is empty.
In Git Bash set `MSYS_NO_PATHCONV=1` for SSM parameter names.
