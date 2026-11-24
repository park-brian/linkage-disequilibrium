# Make an LDlink server copy read 1000G VCFs from /data instead of S3 (the dev bucket has no .tbi indexes).
# Mechanical: "s3://%s/%s" % (aws_info[...bucket...], <path>)  ->  "/data/%s" % (<path>), and checkS3File() is a no-op.
# Everything else is LDlink's code unchanged. usage: python tools/ldlink-local/patch_local_vcf.py <server dir>
import re, sys, pathlib
root = pathlib.Path(sys.argv[1]); total = 0
pat = re.compile(r'"s3://%s/%s"\s*%\s*\(\s*aws_info\[[\'"]bucket[\'"]\]\s*,\s*([^)]+?)\s*\)')
for f in sorted(root.glob("*.py")):
    s = f.open(encoding="utf8", newline="").read(); s2, n = pat.subn(r'"/data/%s" % (\1)', s)  # newline="": keep CRLF/LF as is
    if f.name == "LDcommon.py":
        s2, m = re.subn(r"(def checkS3File\(aws_info, bucket, filePath\):(\r?\n))", r"\1    return True  # local: VCFs are read from /data (patch_local_vcf.py)\2", s2); n += m
    if n: f.open("w", encoding="utf8", newline="").write(s2); total += n; print(f"{f.name}: {n}")
print("total", total)
left = [(f.name, l.strip()) for f in root.glob("*.py") for l in f.read_text(encoding="utf8").splitlines() if "s3://" in l and not l.strip().startswith("#")]
print("remaining s3:// lines:", left)
