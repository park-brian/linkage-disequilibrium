// dbsnp/duplicated_rsids.parquet: rsIDs that have more than one doc in LDlink's dump (296,308). For these, LDlink's
// answer depends on which doc comes first in dump order, so the planner always resolves them through the rsID index.
// usage: node pipeline/dbsnp/dup_rsids.js <releaseRoot>
import { statSync } from "node:fs";
import { zstdCompressSync } from "node:zlib";
import { asyncBufferFromFile, parquetMetadataAsync, parquetReadObjects } from "hyparquet";
import { parquetWriteFile } from "hyparquet-writer";
import { useDataCwd } from "../../tools/config.js";
useDataCwd(); // data paths below are relative to LD_DATA (tools/config.js)
const [root] = process.argv.slice(2);
const file = await asyncBufferFromFile("dbsnp_src/dup_ids.parquet");
const rows = await parquetReadObjects({ file, metadata: await parquetMetadataAsync(file), columns: ["rsid"] });
const ids = BigInt64Array.from(rows, r => BigInt(r.rsid));
for (let i = 1; i < ids.length; i++) if (ids[i] <= ids[i - 1]) throw new Error("not sorted/unique");
const out = `${root}/dbsnp/duplicated_rsids.parquet`;
parquetWriteFile({ filename: out, codec: "ZSTD", compressors: { ZSTD: zstdCompressSync }, rowGroupSize: ids.length, pageSize: 8 << 10,
  columnData: [{ name: "rsid", data: ids, type: "INT64", nullable: false, encoding: "DELTA_BINARY_PACKED", columnIndex: true }] });
console.log(JSON.stringify({ rows: ids.length, bytes: statSync(out).size }));
