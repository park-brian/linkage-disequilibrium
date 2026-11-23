# A local LDlink: the parity oracle

`verify/ldlink/` compares every module port against a real LDlink running on this machine. That LDlink is built
from NCI's own source (`nci-webtools-dceg-linkage`), changed only to read the 1000G VCFs from a local directory.
This folder is the kit; the running copy lives in `$LDLINK_HOME` (default `../ldlink-local`, next to this repo).

**Files:**

| file | what |
|---|---|
| `docker-compose.yml` | `mongo:7` plus LDlink's own backend image (`docker/backend.dockerfile`), on `127.0.0.1:27018` and `127.0.0.1:8090` |
| `patch_local_vcf.py` | the only change to LDlink's code (6 places, 11 lines): `s3://bucket/path` → `/data/path`, and `checkS3File()` returns True |
| `drive.py` | calls LDlink's Python functions directly in the container, for LDhap (warnings the API drops) and LDassoc (no API). The comparers pipe this repo's copy |
| `server.env.example` | LDlink's environment for the container; credential values are placeholders |
| `.env.example` | `MONGO_ROOT_PASSWORD` (a local-only password) and `LDLINK_SRC` (path to the LDlink checkout) |

**Standing it up** (about 6 h, mostly the Mongo restore):

1. **Set up the working copy.** Copy this folder to `$LDLINK_HOME`. Copy `.env.example` → `.env` and
   `server.env.example` → `server.env`, then fill in the placeholders. Every value is local: the app user exists only
   in this container, and the AWS keys can be dummies, because `LDcommon.tabix()` always asks boto3 for credentials.
2. **The server code.** Run `cp -r $LDLINK_SRC/server $LDLINK_HOME/server`, then
   `python -I patch_local_vcf.py $LDLINK_HOME/server`.
3. **The 1000G data.** Under `data/`:
   - `ldlink/1000G/<build>/`: the VCFs and `.tbi` files (hardlinks to `$LD_VCF` are enough);
   - `1000G_population_samples/`, `ldassoc_example/` and `tmp/`;
   - **empty `1000G/<build>/` directories:** LDlink runs tabix with its working directory set to them.
4. **Start Mongo and restore.** Run `docker compose up -d mongo`. Restore the 2023 mongodump with
   `mongorestore --maintainInsertionOrder` (natural order = dump order, which `find_one` parity depends on).
   `gwas_catalog` comes from the 2025-10 backup instead. Skip `api_log`, `*_tmp` and `*_bak`. **Wait for every index
   build to finish**: until then `find_one({id})` is a collection scan. Create the app user from `server.env`.
5. **Start the backend.** Run `docker compose up -d backend`. The mount of `server/` over `/opt/ldlink` hides the
   image's `apache-bin/python3` symlink, so restore it:

   ```
   docker compose exec -u root backend sh -c 'mkdir -p /opt/ldlink/apache-bin && ln -sf /usr/bin/python3.13 /opt/ldlink/apache-bin/python3'
   ```

   `drive.py` puts that directory on `PATH` itself; `docker exec` does not inherit the image's `CMD` environment.
6. **Check it.** `curl "http://127.0.0.1:8090/LDlinkRest/ldpair?var1=rs3&var2=rs4&pop=CEU&genome_build=grch37&json_out=true"`
   should return the rs3/rs4 table. In Git Bash, prefix `docker compose exec` commands with `MSYS_NO_PATHCONV=1`.

The comparers exchange files with the container through `$LDLINK_HOME/data/tmp` (`tools/config.js`: `LDLINK_TMP`).
