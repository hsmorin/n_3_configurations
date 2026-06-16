const STATUS_TEXT = {
    confirmed: 'Confirmed Calabi-Yau',
    confirmed_star: 'Likely Calabi-Yau',
    unconfirmed: 'Unknown Type',
    unconfirmed_star: 'Rational',
    unknown: 'Unknown',
};

function renderMath(tex) {
    try {
        return katex.renderToString(tex, { throwOnError: false, displayMode: false });
    } catch (e) {
        return `<span style="color:red;">${tex}</span>`;
    }
}

function renderText(str) {
    if (!str) return '';
    return str.replace(/\$([^$]+)\$/g, (_, tex) =>
        katex.renderToString(tex, { throwOnError: false, displayMode: false })
    );
}

// Render a small metadata pill: label + value
function metaItem(label, value) {
    if (value === undefined || value === null || value === '') return '';
    // Values containing \ are likely LaTeX
    const rendered = String(value).includes('\\') ? renderMath(value) : value;
    return `
      <div class="meta-item">
        <span class="meta-label">${label}</span>
        <span class="meta-value">${rendered}</span>
      </div>`;
}

const titleEl = document.getElementById('site-title');
titleEl.innerHTML = renderText(titleEl.textContent);

function buildCard(variety) {
    const card = document.createElement('div');
    card.className = 'variety-card';
    card.dataset.id = variety.id;
    card.dataset.status = variety.status ?? 'unknown';

    const statusLabel = STATUS_TEXT[variety.status] ?? '—';
    const hasNotes = variety.notes && variety.notes.trim() !== '';
    const hasHomogenizations = variety.homogenizations &&
        variety.homogenizations.trim() !== '' &&
        variety.homogenizations.trim() !== 'something';

    card.innerHTML = `
    <div class="variety-header">
      <div class="header-left">
        <span class="family-id">${variety.id}</span>
        <span class="variety-name">${renderText(variety.name ?? '')}</span>
      </div>
      <div class="header-right">
        <span class="status-badge">${statusLabel}</span>
        <span class="chevron">&#8964;</span>
      </div>
    </div>
    <div class="variety-body">
      <div class="body-section">

        <div class="meta-grid">
          ${metaItem('Degree', variety.degree)}
          ${metaItem('Dimension', variety.dimension)}
          ${metaItem('Dim. of singular locus', variety['dimension of singular locus'])}
          ${metaItem('Normal in \u{1D538}\u2074', variety['normal in A^4'])}
          ${metaItem('Rational', variety.rational)}
        </div>

        <div class="divider"></div>

        <div class="info-block">
          <p class="block-label">Defining Polynomial</p>
          <div class="math-block">${variety.polynomial && variety.polynomial !== '0'
            ? renderMath(variety.polynomial)
            : '<span class="not-computed">Not yet computed</span>'
        }</div>
        </div>

        ${hasHomogenizations ? `
        <div class="info-block">
          <p class="block-label">Homogenizations</p>
          <div class="math-block">${renderMath(variety.homogenizations)}</div>
        </div>` : ''}

        ${hasNotes ? `
        <div class="divider"></div>
        <div class="info-block">
          <p class="block-label">Notes</p>
          <p class="block-text">${renderText(variety.notes)}</p>
        </div>` : ''}

      </div>
    </div>
  `;

    card.querySelector('.variety-header').addEventListener('click', () => {
        card.classList.toggle('open');
    });

    return card;
}

fetch('data.json')
    .then(res => res.json())
    .then(varieties => {
        varieties.sort((a, b) =>
            a.id.localeCompare(b.id, undefined, { numeric: true })
        );
        const list = document.querySelector('.variety-list');
        varieties.forEach(v => list.appendChild(buildCard(v)));
    })
    .catch(err => console.error('Failed to load data.json:', err));
