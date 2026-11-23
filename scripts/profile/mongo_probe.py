# Read-only probe of the LDLink MongoDB. Run inside the ldlink api container (uses its MONGODB_* env).
# Prints counts, indexes, sample docs, and _id time ranges. Never prints credentials.
import os, json
from pymongo import MongoClient
from bson import ObjectId

e = os.environ
client = MongoClient(host=e["MONGODB_HOST"], port=int(e["MONGODB_PORT"]),
                     username=e["MONGODB_USERNAME"], password=e["MONGODB_PASSWORD"],
                     authSource=e["MONGODB_DATABASE"], serverSelectionTimeoutMS=15000)
db = client[e["MONGODB_DATABASE"]]
out = {"server_version": client.server_info().get("version"), "collections": {}}
skip = {"api_users", "api_log"}  # user PII, not needed
for name in sorted(db.list_collection_names()):
    c = db[name]
    info = {"count": c.estimated_document_count()}
    try:
        info["indexes"] = {k: v["key"] for k, v in c.index_information().items()}
    except Exception as ex:
        info["indexes"] = str(ex)
    if name not in skip:
        info["samples"] = [{k: v for k, v in d.items() if k != "_id"} for d in c.find().limit(2)]
    try:
        first = c.find_one(sort=[("_id", 1)], projection={"_id": 1})
        last = c.find_one(sort=[("_id", -1)], projection={"_id": 1})
        if first and isinstance(first["_id"], ObjectId):
            info["id_time"] = [first["_id"].generation_time.isoformat(), last["_id"].generation_time.isoformat()]
    except Exception:
        pass
    out["collections"][name] = info
print(json.dumps(out, default=str))
