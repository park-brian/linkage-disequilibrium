import time
import threading
from fastparquet import ParquetFile
from guppy import hpy

# filename = "dbsnp100/dbsnp.0.parquet"
filename = "dbsnp2/dbsnp.0.parquet"
# filename = "dbsnp.parquet"
ids = [10468280,10680802,10718688,10852521,11075985,11075987,11075988,11075989,11075990,11075991,11075992,11075993,1121980,113191842,11642015,11642841,11646715,11647020,11649091,12149182,12149574,12149832,1421085,1421086,1477196,1558901,1558902,17817288,17817449,17817497,17817712,17817964,1861866,1861867,2042031,2058908,28429148,28432761,3751812,3751813,3751814,4783819,4783821,4784323,55872725,56094641,56137030,56313538,62033399,62033400,62033403,62033404,62033405,62033406,62033408,62033413,62048402,66908032,7185735,7187250,7187961,7188250,7190396,7190492,7193144,7200298,7201850,7202116,7202296,7206122,7206410,7206629,7206790,72805611,72805612,72805613,79994966,8043757,8044769,8045332,8047395,8047587,8050136,8055197,8057044,8063057,9922047,9922619,9922708,9923147,9923233,9923312,9923544,9926289,9928094,9930333,9930397,9930501,9930506,9931900,9932754,9933040,9933509,9935401,9936385,9937053,9937354,9937709,9939609,9939973,9940128,9940278,9940646,9941349,9972653]
num_partitions = 100
dbs = {}

def get_partition_id(id):
    return id % num_partitions

def load_partition(id):
    # return ParquetFile("dbsnp.parquet")
    return ParquetFile(f"dbsnp2/dbsnp.{get_partition_id(id)}.parquet")

def load_all_partitions():
    for i in range(num_partitions):
        dbs[i] = load_partition(i)

def query_id_by_partition(id):
    partition_id = get_partition_id(id)
    if partition_id not in dbs:
        dbs[partition_id] = load_partition(id)
    partition = dbs[partition_id]
    return partition.to_pandas(['id', 'position_grch37'], filters=[('id', '=', id)], row_filter=True)

def query_ids_by_partition(ids):
    results = []
    for id in ids:
        result = query_id_by_partition(id)
        results.append(result)
    return results


def query_ids_by_partition_threaded(ids):
    threads = []
    for id in ids:
        thread = threading.Thread(target=query_id_by_partition, args=(id,))
        thread.start()
        threads.append(thread)

        # result = query_id_by_partition(id)
        # results.append(result)
    for thread in threads:
        thread.join()

    return threads

filters = [
    ('id', 'in', ids)
    # ('id', '=', 3)
]



# pf = ParquetFile(filename)

start = time.time()
# load_all_partitions()
# h = hpy()
# print(h.heap())

results = query_ids_by_partition(ids)
# df = pf.to_pandas(['id', 'position_grch37'], filters=filters, row_filter=True)
end = time.time()
# print(h.heap())

print("Time: ", end - start)

print(results)