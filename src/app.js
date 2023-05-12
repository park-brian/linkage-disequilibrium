import { getTableOptions, getLinkageDisequilibrium, createTable, exportDelimitedTextFile } from "./utils";
import { exampleSnps } from "./sources";

const loadExampleSnpsButton = document.querySelector("#loadExampleSnpsButton");
loadExampleSnpsButton.addEventListener("click", loadExampleSnps);

const toggleAllPopulationsButton = document.querySelector("#toggleAllPopulationsButton");
toggleAllPopulationsButton.addEventListener("click", toggleAllPopulations);

const resultsDownloadLinks = document.querySelector("#resultsDownloadLinks");
const downloadDprimeResultsButton = document.querySelector("#downloadDprimeResultsButton");
const downloadRsquaredResultsButton = document.querySelector("#downloadRsquaredResultsButton");

const inputForm = document.querySelector("#inputForm");
inputForm.addEventListener("submit", handleSubmit);
inputForm.addEventListener("reset", handleReset);
const resultsElement = document.querySelector("#results");
const resultsTableElement = document.querySelector("#resultsTable");

async function handleSubmit(event) {
  event.preventDefault();
  const form = event.target;
  const start = window.performance.now();

  try {
    form.inputFormFieldset.disabled = true;
    form.submit.disabled = true;
    form.reset.disabled = true;
    resultsElement.innerHTML = "Loading...";
    resultsTableElement.innerHTML = "";
    resultsDownloadLinks.hidden = true;
  
    const genomeBuild = form.genomeBuild.value;
    const snps = form.snps.value.split(/\s+/).filter(Boolean);
    const populations = Array.from(form.populations.selectedOptions).map((option) => option.value);

    const results = await getLinkageDisequilibrium(snps, populations, genomeBuild);
    const tableOptions = getTableOptions(results);

    resultsElement.innerHTML = "";
    resultsTableElement.append(createTable(tableOptions));
    downloadDprimeResultsButton.onclick = () => exportDelimitedTextFile(tableOptions.dPrimeData, "ld-dprime.txt");
    downloadRsquaredResultsButton.onclick = () => exportDelimitedTextFile(tableOptions.rSquaredData, "ld-rsquared.txt");
    resultsDownloadLinks.hidden = false;
  } catch (e) {
    console.error(e);
    resultsElement.innerHTML = e.toString();
  } finally {
    form.inputFormFieldset.disabled = false;
    form.submit.disabled = false;
    form.reset.disabled = false;
    const end = window.performance.now();
    const elapsed = (end - start) / 1000;
    console.log(`Elapsed time: ${elapsed.toFixed(2)} s`);
  }
}

function handleReset(event) {
  resultsElement.innerHTML = "";
  resultsTableElement.innerHTML = "";
  resultsDownloadLinks.hidden = true;
}

function loadExampleSnps() {
  const snpsElement = document.querySelector("#snps");
  snpsElement.value = exampleSnps.join("\n");
}

function toggleAllPopulations() {
  const populationsElement = document.querySelector("#populations");
  const options = Array.from(populationsElement.options);
  const allSelected = options.every((option) => option.selected);
  options.forEach((option) => option.selected = !allSelected);
}