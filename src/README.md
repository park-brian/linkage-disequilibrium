# src/: the browser-ready core

**What it is:** runtime code for browser and Node, with no build step.
- **Dependencies:** `hyparquet` and `hyparquet-compressors` only.
- **Data access:** a caller-supplied `open(path)` reads the release ([docs/DATA.md](../docs/DATA.md)).

**More detail:**
- the module map and the one-way dependency rule: [docs/ARCHITECTURE.md](../docs/ARCHITECTURE.md);
- what the LDlink ports reproduce: [docs/PARITY.md](../docs/PARITY.md);
- the rules for changing anything here: [AGENTS.md](../AGENTS.md).

**Using a module port:**

```js
import { loadKernel, KERNEL_WASM } from "./ld/kernel.mjs";
import { createContext } from "./ldlink/context.mjs";
import { ldpair } from "./ldlink/ldpair.mjs";
const kernel = await loadKernel(new Uint8Array(await (await fetch(KERNEL_WASM)).arrayBuffer()));
const ctx = await createContext({ root, open, readPanel, kernel }); // readPanel(pop) → sample IDs
const out = await ldpair(ctx, { var1: "rs3", var2: "rs4", pop: "CEU", build: "grch37" });
```
