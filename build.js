import { build, context } from "esbuild";
import { polyfillNode } from "esbuild-plugin-polyfill-node";
const args = process.argv.slice(2);

const config = {
  entryPoints: ["src/app.js"],
  bundle: true,
  minify: true,
  sourcemap: true,
  outfile: "docs/app.js",
  plugins: [
    polyfillNode({
      polyfills: ["buffer"],
    }),
  ],
}

if (args.includes("watch")) {
  const ctx = await context(config);
  await ctx.watch();
  const { port } = await ctx.serve({
    servedir: 'docs',
  })
  console.log(`Serving on http://localhost:${port}`);
} else {
  build(config);
}
