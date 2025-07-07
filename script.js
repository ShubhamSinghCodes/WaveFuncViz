document.addEventListener('DOMContentLoaded', function () {
    const plotButton = document.getElementById('plotButton');
    const nInput = document.getElementById('n');
    const lInput = document.getElementById('l');
    const mInput = document.getElementById('m');
    const orbitalNameEl = document.getElementById('orbitalName');
    const wavefunctionEl = document.getElementById('wavefunction');
    const plotEl = document.getElementById('plot');
    const viewToggle = document.getElementById('viewToggle');

    let useIsosurface = false;

    // --- Custom Simplification Rules for math.js ---
    const customRules = [
        // Move constants to the left of variables
        { l: 'v * c', r: 'c * v'},
        // Flatten nested divisions
        { l: '(n1 / n2) / n3', r: 'n1 / (n2 * n3)'},
        // Explicitly define Pythagorean identities
        { l: '1 - cos(n)^2', r: 'sin(n)^2'},
        { l: '1 - sin(n)^2', r: 'cos(n)^2'}
    ].concat(math.simplify.rules);
    const symbolicScope = {
        a_0: 1,
        pi: 3.14
    };


    // --- Mathematical Functions (Numeric for Plotting) ---
    const factorial = (n) => {
        if (n < 0 || n % 1 !== 0) return NaN;
        if (n === 0) return 1;
        let res = 1;
        for (let i = 2; i <= n; i++) res *= i;
        return res;
    };

    const binomialCoefficient = (n, k) => {
        if (k < 0 || k > n) return 0;
        return factorial(n) / (factorial(k) * factorial(n - k));
    };

    const associatedLaguerre = (p, q, x) => {
        let sum = 0;
        for (let i = 0; i <= p; i++) {
            sum += Math.pow(-1, i) * binomialCoefficient(p + q, p - i) * Math.pow(x, i) / factorial(i);
        }
        return sum;
    };

    const stableAssociatedLegendre = (l, m, x) => {
        m = Math.abs(m);
        let pmm = 1.0;
        if (m > 0) {
            const sqrt_one_minus_x2 = Math.sqrt(1.0 - x * x);
            for (let i = 1; i <= m; i++) {
                pmm *= (-1.0) * (2 * i - 1) * sqrt_one_minus_x2;
            }
        }
        if (l === m) return pmm;
        let pmm1 = x * (2 * m + 1) * pmm;
        if (l === m + 1) return pmm1;
        let pll = 0;
        for (let i = m + 2; i <= l; i++) {
            pll = (x * (2 * i - 1) * pmm1 - (i + m - 1) * pmm) / (i - m);
            pmm = pmm1;
            pmm1 = pll;
        }
        return pll;
    };

    const sphericalHarmonic = (m, l, theta, phi) => {
        const abs_m = Math.abs(m);
        let normalization = Math.sqrt(((2 * l + 1) / (4 * Math.PI)) * (factorial(l - abs_m) / factorial(l + abs_m)));
        if (m !== 0) {
            normalization *= Math.sqrt(2);
        }
        const p = stableAssociatedLegendre(l, abs_m, Math.cos(theta));
        if (m >= 0) {
            return normalization * p * Math.cos(m * phi);
        } else {
            return normalization * p * Math.sin(abs_m * phi);
        }
    };

    const radialWavefunction = (n, l, r) => {
        const a0 = 1;
        const rho = 2 * r / (n * a0);
        const normalization = Math.sqrt(Math.pow(2 / (n * a0), 3) * (factorial(n - l - 1) / (2 * n * factorial(n + l))));
        const laguerre = associatedLaguerre(n - l - 1, 2 * l + 1, rho);
        return normalization * Math.exp(-rho / 2) * Math.pow(rho, l) * laguerre;
    };

    const wavefunction = (n, l, m, r, theta, phi) => {
        const R = radialWavefunction(n, l, r);
        if (isNaN(R)) return 0;
        const Y = sphericalHarmonic(m, l, theta, phi);
        return R * Y;
    };

    // --- UI and Plotting ---

    const getOrbitalName = (n, l) => {
        const subshell = ['s', 'p', 'd', 'f', 'g', 'h', 'i'][l] || `l=${l}`;
        return `${n}${subshell}`;
    };

    // --- Symbolic Polynomial Generation (as math.js strings) ---
    const symbolicLaguerreMathJS = (p, q, varStr) => {
        if (p === 0) return "1";
        let sum_str = "";
        for (let i = 0; i <= p; i++) {
            const coeff_str = `((-1)^${i} * combinations(${p + q}, ${p - i})) / (${i}!)`;
            const term_str = `(${coeff_str}) * ${varStr}^${i}`;
            sum_str += (i > 0 ? " + " : "") + term_str;
        }
        return `(${sum_str})`;
    };

    const symbolicLegendreMathJS = (l, m, varStr) => {
        const abs_m = Math.abs(m);
        let sum_str = "";
        const floor_val = Math.floor((l - abs_m) / 2);

        for (let k = 0; k <= floor_val; k++) {
            const coeff_str = `((-1)^${k} * (2*${l} - 2*${k})!) / (${k}! * (${l}-${k})! * (${l} - ${abs_m} - 2*${k})!)`;
            const power = l - abs_m - 2 * k;
            const term_str = `(${coeff_str}) * ${varStr}^${power}`;
            sum_str += (k > 0 ? " + " : "") + term_str;
        }

        const prefactor_str = `(1 / 2^${l}) * (1 - ${varStr}^2)^(${abs_m}/2)`;
        return `${prefactor_str} * (${sum_str})`;
    };

    const generateWavefunctionEquation = (n, l, m) => {
        const abs_m = Math.abs(m);
        const rho_def = `(2*r/(${n}*a_0))`;

        try {
            // --- Build all parts of the expression as strings ---
            const N_nl = `sqrt((2/(${n}*a_0))^3 * (${n-l-1})! / (2*${n}*(${n+l})!))`;

            const laguerre_expr_str = symbolicLaguerreMathJS(n - l - 1, 2 * l + 1, rho_def);

            const R_part = `${N_nl} * e^(-${rho_def}/2) * ${rho_def}^${l} * (${laguerre_expr_str})`;

            const N_lm = `sqrt((${m === 0 ? '1' : '2'} * (2*${l}+1)) / (4*pi) * ((${l-abs_m})!) / ((${l+abs_m})!))`;
            const legendre_expr_str = symbolicLegendreMathJS(l, m, "(cos(theta))");
            const phi_part = m >= 0 ? `cos(${m}*phi)` : `sin(${abs_m}*phi)`;
            const Y_part = `${N_lm} * ${legendre_expr_str} * ${phi_part}`;

            // --- Combine everything and simplify with the custom ruleset ---
            const full_expr_str = `(${R_part}) * (${Y_part})`;
            const simplified_node = math.simplify(full_expr_str, customRules, symbolicScope);

            const final_latex_expr = simplified_node.toTex();

            katex.render(`\\psi_{${n},${l},${m}} = ${final_latex_expr}`, wavefunctionEl, {
                displayMode: true,
                throwOnError: false
            });

        } catch(e) {
            wavefunctionEl.textContent = "Error processing equation with math.js.";
            console.error(e);
        }
    };

    const updateUI = () => {
        const n = parseInt(nInput.value);
        const l = parseInt(lInput.value);
        const m = parseInt(mInput.value);

        if (n <= 0 || l < 0 || l >= n || Math.abs(m) > l) {
            alert("Invalid quantum numbers! Please ensure n > 0, 0 ≤ l < n, and -l ≤ m ≤ l.");
            return;
        }

        orbitalNameEl.innerHTML = getOrbitalName(n, l) + " (In the units where a<sub>0</sub> is 1)";
        generateWavefunctionEquation(n, l, m);
        plotOrbital(n, l, m);
    };

    const plotOrbital = (n, l, m) => {
        const gridSize = 40;
        const range = n > 1 ? n * n * 2.75 : 5;

        const x_full = [], y_full = [], z_full = [], wave_values_full = [];

        for (let i = 0; i < gridSize; i++) {
            for (let j = 0; j < gridSize; j++) {
                for (let k = 0; k < gridSize; k++) {
                    const xi = -range + 2 * range * i / (gridSize - 1);
                    const yi = -range + 2 * range * j / (gridSize - 1);
                    const zi = -range + 2 * range * k / (gridSize - 1);

                    const r = Math.sqrt(xi * xi + yi * yi + zi * zi);
                    if (r === 0) continue;

                    const theta = Math.acos(zi / r);
                    const phi = Math.atan2(yi, xi);

                    const psi = wavefunction(n, l, m, r, theta, phi);

                    x_full.push(xi);
                    y_full.push(yi);
                    z_full.push(zi);
                    wave_values_full.push(psi);
                }
            }
        }

        const traces = [];

        if (useIsosurface) {
            const max_abs_val = Math.max(...wave_values_full.map(v => Math.abs(v)));
            const threshold = max_abs_val * 0.1;
            traces.push({
                type: 'isosurface',
                x: x_full, y: y_full, z: z_full,
                value: wave_values_full,
                opacity: 1,
                isomin: threshold,
                isomax: max_abs_val,
                colorscale: [[0, 'rgb(0,0,255)'], [1, 'rgb(0,0,255)']],
                showscale: false,
                name: 'Positive Lobe'
            });
            traces.push({
                type: 'isosurface',
                x: x_full, y: y_full, z: z_full,
                value: wave_values_full,
                opacity: 1,
                isomin: -max_abs_val,
                isomax: -threshold,
                colorscale: [[0, 'rgb(255,0,0)'], [1, 'rgb(255,0,0)']],
                showscale: false,
                name: 'Negative Lobe'
            });
        } else {
            const wave_prob_full = wave_values_full.map(v => v*v);
            const max_prob = Math.max(...wave_prob_full);
            const threshold = 0.01 * max_prob;

            const x = [], y = [], z = [], wave_values = [];

            for(let i=0; i < wave_prob_full.length; i++) {
                if(wave_prob_full[i] > threshold) {
                    x.push(x_full[i]);
                    y.push(y_full[i]);
                    z.push(z_full[i]);
                    wave_values.push(wave_values_full[i]);
                }
            }

            const wave_prob = wave_values.map(v => v*v);
            const current_max_prob = Math.max(...wave_prob);
            const sizes = wave_prob.map(p => Math.min(5, 150 * (p / current_max_prob)));
            const colors = wave_values.map(v => v >= 0 ? 'blue' : 'red');

            traces.push({
                x, y, z,
                mode: 'markers',
                type: 'scatter3d',
                marker: {
                    color: colors,
                    size: sizes,
                    opacity: 1,
                    sizemin: 2,
                },
            });
        }

        const layout = {
            title: `${getOrbitalName(n, l)} Orbital`,
            margin: { l: 0, r: 0, b: 20, t: 40 },
             scene: {
                xaxis: { title: 'x', showgrid: false, zeroline: false, showline: false, ticks: '', showticklabels: false },
                yaxis: { title: 'y', showgrid: false, zeroline: false, showline: false, ticks: '', showticklabels: false },
                zaxis: { title: 'z', showgrid: false, zeroline: false, showline: false, ticks: '', showticklabels: false },
                camera: { eye: { x: 1.5, y: 1.5, z: 1.5 } },
                aspectratio: { x: 1, y: 1, z: 1 },
            },
            paper_bgcolor: 'rgba(0,0,0,0)',
            plot_bgcolor: 'rgba(0,0,0,0)'
        };

        Plotly.newPlot(plotEl, traces, layout, {responsive: true});
    };

    // Event Listeners
    plotButton.addEventListener('click', updateUI);
    viewToggle.addEventListener('change', (e) => {
        useIsosurface = e.target.checked;
        updateUI();
    });

    // Initial plot
    updateUI();
});