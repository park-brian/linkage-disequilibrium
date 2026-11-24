# src/: the browser-ready core

**What it is:** runtime code for browser and Node, with no build step.
- **Dependencies:** `hyparquet` and `hyparquet-compressors` only.
- **Data access:** a caller-supplied `open(path)` reads the release ([docs/DATA.md](../docs/DATA.md)).

**More detail:**
- the module map and the one-way dependency rule: [docs/ARCHITECTURE.md](../docs/ARCHITECTURE.md);
- what the LDlink ports reproduce: [docs/PARITY.md](../docs/PARITY.md);
- the rules for changing anything here: [AGENTS.md](../AGENTS.md).

**The LD statistics** (the core: numbers and typed errors, allele-aware rows; `core/ld.js`):

```js
import { createLD } from "./core/ld.js";
const ld = createLD({ root, open, kernel });
const { r2, dprime, haplotypes } = await ld.pair("rs3", "rs4", { pop: "CEU" }); // build: "grch37" by default
const { proxies } = await ld.proxy("rs3", { pop: "CEU", window: 500000, minR2: 0.8 });
```

**LDlink's exact output** (a module port: `parity/ldlink/`):

```js
import { loadKernel, KERNEL_WASM } from "./core/ld/kernel.js";
import { createContext } from "./parity/ldlink/context.js";
import { ldpair } from "./parity/ldlink/ldpair.js";
const kernel = await loadKernel(new Uint8Array(await (await fetch(KERNEL_WASM)).arrayBuffer()));
const ctx = await createContext({ root, open, readPanel, kernel }); // readPanel(pop) → sample IDs
const out = await ldpair(ctx, { var1: "rs3", var2: "rs4", pop: "CEU", build: "grch37" });
```
