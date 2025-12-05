// HGT-Detector Frontend Logic

// Configuration
const CONFIG = {
    apiBaseUrl: '/api/v1'
};

// DOM Elements
const tabs = document.querySelectorAll('.tab-btn');
const tabPanes = document.querySelectorAll('.tab-pane');
const dropZone = document.getElementById('drop-zone');
const fileInput = document.getElementById('file-input');
const fileNameDisplay = document.getElementById('file-name');
const modelSelect = document.getElementById('model-select');
const modelDesc = document.getElementById('model-description');
const analyzeBtn = document.getElementById('analyze-btn');
const dashboard = document.getElementById('dashboard');

// State
let selectedFile = null;

// --- Event Listeners ---

// Tab Switching
tabs.forEach(tab => {
    tab.addEventListener('click', () => {
        tabs.forEach(t => t.classList.remove('active'));
        tabPanes.forEach(p => p.classList.remove('active'));

        tab.classList.add('active');
        document.getElementById(tab.dataset.tab).classList.add('active');
    });
});

// File Upload (Drag & Drop)
dropZone.addEventListener('dragover', (e) => {
    e.preventDefault();
    dropZone.classList.add('dragover');
});

dropZone.addEventListener('dragleave', () => {
    dropZone.classList.remove('dragover');
});

dropZone.addEventListener('drop', (e) => {
    e.preventDefault();
    dropZone.classList.remove('dragover');
    if (e.dataTransfer.files.length) {
        handleFile(e.dataTransfer.files[0]);
    }
});

fileInput.addEventListener('change', (e) => {
    if (e.target.files.length) {
        handleFile(e.target.files[0]);
    }
});

function handleFile(file) {
    if (file.name.endsWith('.fasta') || file.name.endsWith('.fa')) {
        selectedFile = file;
        fileNameDisplay.textContent = `Selected: ${file.name}`;
    } else {
        alert('Please upload a valid .fasta or .fa file.');
    }
}

// Model Selection Description Update
const modelDescriptions = {
    near: "Detects transfer between closely related species (e.g., E. coli to Legionella).",
    moderate: "Detects transfer between moderately distant species (e.g., Soil Bacteria to Host).",
    distant: "Detects transfer between distantly related species (e.g., Eukaryote to Bacteria)."
};

modelSelect.addEventListener('change', (e) => {
    modelDesc.textContent = modelDescriptions[e.target.value];
});

// Analyze Button
analyzeBtn.addEventListener('click', runAnalysis);

// --- Analysis Logic ---

async function runAnalysis() {
    const activeTab = document.querySelector('.tab-btn.active').dataset.tab;

    analyzeBtn.textContent = 'Analyzing...';
    analyzeBtn.disabled = true;

    try {
        const formData = new FormData();
        formData.append('distance', modelSelect.value);
        formData.append('foreign_threshold', '0.05');

        if (activeTab === 'upload-tab') {
            if (!selectedFile) {
                alert('Please select a file first.');
                analyzeBtn.textContent = 'Run Analysis';
                analyzeBtn.disabled = false;
                return;
            }
            formData.append('file', selectedFile);
        } else {
            const fastaText = document.getElementById('sequence-input').value;
            if (!fastaText.trim()) {
                alert('Please paste a sequence.');
                analyzeBtn.textContent = 'Run Analysis';
                analyzeBtn.disabled = false;
                return;
            }
            formData.append('fasta', fastaText);
        }

        const resp = await fetch(`${CONFIG.apiBaseUrl}/predict`, {
            method: 'POST',
            body: formData
        });
        const data = await resp.json();
        if (!resp.ok) {
            throw new Error(data.error || 'Prediction failed');
        }
        displayResults(data);
        dashboard.classList.remove('hidden');
        dashboard.scrollIntoView({ behavior: 'smooth' });
    } catch (err) {
        alert(err.message || 'Error during analysis');
        console.error(err);
    } finally {
        analyzeBtn.textContent = 'Run Analysis';
        analyzeBtn.disabled = false;
    }
}

function displayResults(data) {
    // Display counts
    document.getElementById('result-total').textContent = data.totalGenes;
    document.getElementById('result-host').textContent = `${data.hostCount} (${data.hostPercent}%)`;
    document.getElementById('result-foreign').textContent = `${data.foreignCount} (${data.foreignPercent}%)`;

    // Render Charts
    renderPCA(data.pcaData);
    renderHeatmap(data.heatmapData);
    renderBarChart(data.barData);
}

// --- Chart Rendering (Chart.js) ---

let pcaChartInstance = null;
let heatmapChartInstance = null;
let barChartInstance = null;

function renderPCA(data) {
    const ctx = document.getElementById('pcaChart').getContext('2d');

    if (pcaChartInstance) pcaChartInstance.destroy();

    pcaChartInstance = new Chart(ctx, {
        type: 'scatter',
        data: {
            datasets: [
                {
                    label: 'Host Reference',
                    data: data.host,
                    backgroundColor: 'rgba(0, 128, 128, 0.5)',
                    pointRadius: 4
                },
                {
                    label: 'Donor Reference',
                    data: data.donor,
                    backgroundColor: 'rgba(231, 76, 60, 0.5)',
                    pointRadius: 4
                },
                {
                    label: 'Input (Pred. Host)',
                    data: data.inputHost || [],
                    backgroundColor: '#2ecc71',
                    pointRadius: 6,
                    pointStyle: 'triangle'
                },
                {
                    label: 'Input (Pred. Foreign)',
                    data: data.inputForeign || [],
                    backgroundColor: '#9b59b6',
                    pointRadius: 6,
                    pointStyle: 'triangle'
                }
            ]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            scales: {
                x: { title: { display: true, text: 'PC1 (Codon Usage)' } },
                y: { title: { display: true, text: 'PC2 (Codon Usage)' } }
            }
        }
    });
}

function renderHeatmap(data) {
    const ctx = document.getElementById('heatmapChart').getContext('2d');

    if (heatmapChartInstance) heatmapChartInstance.destroy();

    heatmapChartInstance = new Chart(ctx, {
        type: 'line',
        data: {
            labels: data.map((_, i) => i),
            datasets: [{
                label: 'Foreign Probability',
                data: data,
                borderColor: '#e74c3c',
                backgroundColor: 'rgba(231, 76, 60, 0.2)',
                fill: true,
                tension: 0.4
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            scales: {
                y: {
                    beginAtZero: true,
                    max: 1,
                    title: { display: true, text: 'Probability (Foreign)' }
                },
                x: {
                    title: { display: true, text: 'Gene Index' }
                }
            },
            plugins: {
                legend: { display: false }
            }
        }
    });
}

function renderBarChart(data) {
    const ctx = document.getElementById('barChart').getContext('2d');

    if (barChartInstance) barChartInstance.destroy();

    barChartInstance = new Chart(ctx, {
        type: 'bar',
        data: {
            labels: data.map(d => d.gene),
            datasets: [{
                label: 'Foreign Probability',
                data: data.map(d => d.prob),
                backgroundColor: data.map(d => d.prob > 0.5 ? '#e74c3c' : '#008080')
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            scales: {
                y: {
                    beginAtZero: true,
                    max: 1,
                    title: { display: true, text: 'Probability' }
                }
            }
        }
    });
}
