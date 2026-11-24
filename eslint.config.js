// Correctness only: names that don't exist, values never used, code that can't run, and the like. No style rules: the
// code's density is deliberate (AGENTS.md). Run by npm test (test/lint.test.js); npx eslint . lists the findings.
import globals from "globals";

const correctness = {
  "no-undef": "error",
  "no-unused-vars": ["error", { args: "after-used", caughtErrors: "none", varsIgnorePattern: "^_", argsIgnorePattern: "^_",
    ignoreRestSiblings: true }],
  "no-unreachable": "error",
  "no-dupe-keys": "error",
  "no-dupe-args": "error",
  "no-dupe-else-if": "error",
  "no-duplicate-case": "error",
  "no-self-assign": "error",
  "no-self-compare": "error",
  "no-constant-condition": ["error", { checkLoops: false }],
  "no-constant-binary-expression": "error",
  "no-cond-assign": ["error", "except-parens"],
  "no-const-assign": "error",
  "no-func-assign": "error",
  "no-import-assign": "error",
  "no-redeclare": "error",
  "no-shadow-restricted-names": "error",
  "no-loss-of-precision": "error",
  "no-unsafe-finally": "error",
  "no-unsafe-negation": "error",
  "no-unsafe-optional-chaining": "error",
  "no-sparse-arrays": "error",
  "no-empty-pattern": "error",
  "no-obj-calls": "error",
  "no-compare-neg-zero": "error",
  "use-isnan": "error",
  "valid-typeof": "error",
  "for-direction": "error",
  "getter-return": "error",
};

export default [
  { ignores: ["node_modules/**", "dist/**", "docs/history/**"] },
  { files: ["**/*.js"], languageOptions: { ecmaVersion: "latest", sourceType: "module" }, rules: correctness },
  // the runtime: browser and Node alike, so only what both have
  { files: ["src/**/*.js"], languageOptions: { globals: { ...globals["shared-node-browser"] } } },
  { files: ["app/**/*.js", "bench/app_requests.js"], languageOptions: { globals: { ...globals.browser } } },
  { files: ["app/workers/**/*.js"], languageOptions: { globals: { ...globals.worker } } },
  { files: ["pipeline/**/*.js", "verify/**/*.js", "tools/**/*.js", "test/**/*.js", "bench/**/*.js", "eslint.config.js"],
    ignores: ["bench/app_requests.js"], languageOptions: { globals: { ...globals.node } } },
  // their page.evaluate callbacks run in the browser
  { files: ["verify/app/run.js", "verify/app/set.js"], languageOptions: { globals: { ...globals.node, ...globals.browser } } },
  // W. J. Cody's erfc coefficients, digit for digit as published (more digits than a double holds; JS parses them to
  // the same doubles the reference does), so the source can be checked against the paper
  { files: ["src/core/stats.js"], rules: { "no-loss-of-precision": "off" } },
];
