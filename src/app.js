import { getTableOptions, getLinkageDisequilibrium, createTable } from "./utils";
import { exampleSnps } from "./sources";

const loadExampleSnpsButton = document.querySelector("#loadExampleSnpsButton");
loadExampleSnpsButton.addEventListener("click", loadExampleSnps);

const toggleAllPopulationsButton = document.querySelector("#toggleAllPopulationsButton");
toggleAllPopulationsButton.addEventListener("click", toggleAllPopulations);

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
    resultsElement.innerHTML = "";
    resultsTableElement.innerHTML = "";
  
    const genomeBuild = form.genomeBuild.value;
    const snps = form.snps.value.split(/\s+/).filter(Boolean);
    const populations = Array.from(form.populations.selectedOptions).map((option) => option.value);

    const results = await getLinkageDisequilibrium(snps, populations, genomeBuild);
    resultsTableElement.innerHTML = "";
    resultsTableElement.append(createTable(getTableOptions(results)));
    console.log(results);
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
  console.log("reset");
  resultsElement.innerHTML = "";
  resultsTableElement.innerHTML = "";
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