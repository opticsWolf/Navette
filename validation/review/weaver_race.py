# Step 4: OpticalWeaver LRU-plan eviction + thread race.
import threading
import numpy as np
from navette.spectralweave import OpticalWeaver

OK = True
def report(name, cond, detail=""):
    global OK
    print(f"  {name}: {'OK' if cond else 'FAIL'} {detail}")
    OK = bool(cond) and OK

FULL = np.arange(0.0, 600.0, 1.0)
TILES = [(0, 150), (150, 300), (300, 450), (450, 600)]   # 4 tiling frames


def make_weaver(cache_size):
    wv = OpticalWeaver(cache_size=cache_size)
    seed = (0.0, "seed", "x")
    for lo, hi in TILES:
        b = FULL[lo:hi]
        wv.set_data(seed, np.zeros(hi - lo), b)
    return wv


print("=== LRU plan eviction (cache_size=2, 4 grids, 10 keys) ===")
wv = make_weaver(2)
rng = np.random.default_rng(11)
keys = []
for k in range(10):
    key = (550.0 + k, "R", "s")
    curve = rng.random(FULL.size)
    n = wv.unweave(key, FULL, curve)
    if n != len(TILES):
        report("unweave touches all frames", False, f"k={k} n={n}")
    keys.append((key, curve))
bad = 0
for key, curve in keys:
    wl_out, v_out = wv.get_weaved(key)
    if wl_out.size != FULL.size or not np.allclose(v_out, curve, atol=1e-12):
        bad += 1
report("reassembly exact after plan-thrash", bad == 0, f"bad={bad}/10")

print("=== thread race: 8 threads, mixed ops, one weaver ===")
wv2 = make_weaver(2)   # tiny cache on purpose -> maximal plan thrash
errors = []
written = {}   # key -> curve (single writer per key)
wlock = threading.Lock()
def worker(tid):
    try:
        for rep in range(40):
            key = (700.0 + tid * 100 + rep, "T", "p")   # unique per (tid, rep)
            op = rep % 5
            if op in (0, 1):
                curve = np.sin(FULL / (7.0 + tid)) + tid
                wv2.unweave(key, FULL, curve)
                with wlock:
                    written[key] = curve
            elif op == 2:
                curve = np.cos(FULL / (3.0 + tid))
                wv2.unweave_collection(FULL, {key: curve})
                with wlock:
                    written[key] = curve
            elif op == 3:
                try:
                    wv2.get_weaved(key)
                except ValueError:
                    pass   # not written yet by THIS thread's earlier ops
            else:
                wv2.get_weaved_collections()
                if rep % 13 == 0:
                    wv2.invalidate_cache()
    except Exception as e:
        errors.append(f"t{tid}: {type(e).__name__}: {e}")

threads = [threading.Thread(target=worker, args=(t,)) for t in range(8)]
for t in threads: t.start()
for t in threads: t.join(timeout=120)
report("no exceptions under 8-thread race", not errors, str(errors[:2]))

# integrity: every key must reassemble EXACTLY to its single writer's curve
bad = 0
for key, curve in written.items():
    try:
        _, v_out = wv2.get_weaved(key)
        if not np.allclose(v_out, curve, atol=1e-12):
            bad += 1
    except Exception:
        bad += 1
report(f"post-race integrity ({len(written)} single-writer keys)", bad == 0,
       f"torn/mismatched={bad}")

print()
print("ALL OK" if OK else "FAILURES PRESENT")
