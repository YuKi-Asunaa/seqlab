/* ============================================================
   SeqLab · 序列分析工作台
   纯前端实现：序列统计 / 六框 ORF / 限制酶图谱 / 密码子偏好 / 双序列比对
   ============================================================ */
'use strict';

/* ---------------- 基础工具 ---------------- */
const $ = (s, r = document) => r.querySelector(s);
const $$ = (s, r = document) => Array.from(r.querySelectorAll(s));
const clamp = (v, a, b) => (v < a ? a : v > b ? b : v);
const fmt = (n) => Number(n).toLocaleString('en-US');
const pctS = (v, d = 1) => (v * 100).toFixed(d) + '%';

function canvasCtx(canvas) {
  const rect = canvas.parentElement.getBoundingClientRect();
  const w = Math.max(1, Math.round(rect.width));
  const h = Math.max(1, Math.round(rect.height));
  const dpr = Math.min(2, window.devicePixelRatio || 1);
  canvas.width = Math.round(w * dpr);
  canvas.height = Math.round(h * dpr);
  canvas.style.width = w + 'px';
  canvas.style.height = h + 'px';
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
  return { ctx, w, h };
}

function niceStep(total, target) {
  const raw = Math.max(1, total) / Math.max(1, target);
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const norm = raw / mag;
  let step = mag;
  if (norm > 5) step = 10 * mag;
  else if (norm > 2) step = 5 * mag;
  else if (norm > 1) step = 2 * mag;
  return Math.max(1, Math.round(step));
}

function mulberry32(a) {
  return function () {
    a |= 0; a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/* 提示气泡 */
const tip = $('#tooltip');
function showTip(html, x, y) {
  tip.innerHTML = html;
  tip.classList.add('show');
  const r = tip.getBoundingClientRect();
  let left = x + 14, top = y - r.height - 10;
  if (left + r.width > window.innerWidth - 8) left = x - r.width - 14;
  if (top < 8) top = y + 18;
  tip.style.left = left + 'px';
  tip.style.top = top + 'px';
}
function hideTip() { tip.classList.remove('show'); }

/* ---------------- 遗传密码 ---------------- */
const BASE_ORDER = 'TCAG';
const CODONS = [];
for (const b1 of BASE_ORDER) for (const b2 of BASE_ORDER) for (const b3 of BASE_ORDER) CODONS.push(b1 + b2 + b3);
const AA_TABLE = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
const CODON_AA = Object.create(null);
CODONS.forEach((c, i) => { CODON_AA[c] = AA_TABLE[i]; });

const AA_NAMES = {
  A: 'Ala', R: 'Arg', N: 'Asn', D: 'Asp', C: 'Cys', Q: 'Gln', E: 'Glu', G: 'Gly',
  H: 'His', I: 'Ile', L: 'Leu', K: 'Lys', M: 'Met', F: 'Phe', P: 'Pro', S: 'Ser',
  T: 'Thr', W: 'Trp', Y: 'Tyr', V: 'Val', '*': 'Stop',
};
const SYN = Object.create(null);   // aa -> [codon...]
for (const c of CODONS) {
  const aa = CODON_AA[c];
  if (!SYN[aa]) SYN[aa] = [];
  SYN[aa].push(c);
}
const AA20 = Object.keys(SYN).filter((a) => a !== '*').sort();

/* ---------------- 限制性内切酶 ---------------- */
const ENZYMES = [
  { n: 'EcoRI', s: 'GAATTC', cut: 1, c: 1 },
  { n: 'BamHI', s: 'GGATCC', cut: 1, c: 1 },
  { n: 'HindIII', s: 'AAGCTT', cut: 1, c: 1 },
  { n: 'XhoI', s: 'CTCGAG', cut: 1, c: 1 },
  { n: 'SalI', s: 'GTCGAC', cut: 1, c: 1 },
  { n: 'PstI', s: 'CTGCAG', cut: 5, c: 1 },
  { n: 'KpnI', s: 'GGTACC', cut: 5, c: 1 },
  { n: 'SacI', s: 'GAGCTC', cut: 5, c: 1 },
  { n: 'XbaI', s: 'TCTAGA', cut: 1, c: 1 },
  { n: 'SphI', s: 'GCATGC', cut: 5, c: 1 },
  { n: 'SmaI', s: 'CCCGGG', cut: 3, c: 1 },
  { n: 'EcoRV', s: 'GATATC', cut: 3, c: 1 },
  { n: 'PvuII', s: 'CAGCTG', cut: 3, c: 1 },
  { n: 'NheI', s: 'GCTAGC', cut: 1, c: 1 },
  { n: 'SpeI', s: 'ACTAGT', cut: 1, c: 1 },
  { n: 'BglII', s: 'AGATCT', cut: 1, c: 1 },
  { n: 'NcoI', s: 'CCATGG', cut: 1, c: 1 },
  { n: 'NdeI', s: 'CATATG', cut: 2, c: 1 },
  { n: 'NotI', s: 'GCGGCCGC', cut: 2, c: 1 },
  { n: 'MluI', s: 'ACGCGT', cut: 1, c: 1 },
  { n: 'AgeI', s: 'ACCGGT', cut: 1, c: 1 },
  { n: 'ClaI', s: 'ATCGAT', cut: 2, c: 1 },
  { n: 'AatII', s: 'GACGTC', cut: 5, c: 1 },
  { n: 'AvrII', s: 'CCTAGG', cut: 1, c: 1 },
  { n: 'BstBI', s: 'TTCGAA', cut: 3, c: 0 },
  { n: 'AflII', s: 'CTTAAG', cut: 5, c: 0 },
  { n: 'BspEI', s: 'TCCGGA', cut: 3, c: 0 },
  { n: 'NsiI', s: 'ATGCAT', cut: 5, c: 0 },
  { n: 'ScaI', s: 'AGTACT', cut: 3, c: 0 },
  { n: 'StuI', s: 'AGGCCT', cut: 3, c: 0 },
  { n: 'PvuI', s: 'CGATCG', cut: 4, c: 0 },
  { n: 'BspHI', s: 'TCATGA', cut: 1, c: 0 },
  { n: 'AluI', s: 'AGCT', cut: 2, c: 0, four: 1 },
  { n: 'HaeIII', s: 'GGCC', cut: 2, c: 0, four: 1 },
  { n: 'MspI', s: 'CCGG', cut: 1, c: 0, four: 1 },
  { n: 'TaqI', s: 'TCGA', cut: 1, c: 0, four: 1 },
  { n: 'Sau3AI', s: 'GATC', cut: 0, c: 0, four: 1 },
  { n: 'RsaI', s: 'GTAC', cut: 2, c: 0, four: 1 },
  { n: 'DpnII', s: 'GATC', cut: 0, c: 0, four: 1 },
];

function enzymeEnd(e) {
  const L = e.s.length;
  if (e.cut === L / 2) return { t: '平末端', cls: 'four' };
  return e.cut < L / 2 ? { t: '5′突出', cls: 'six' } : { t: '3′突出', cls: 'six' };
}

function findSites(seq, enzymeList, circular) {
  const sites = [];
  const len = seq.length;
  for (const e of enzymeList) {
    const L = e.s.length;
    if (len < L) continue;
    const hay = circular ? seq + seq.slice(0, L - 1) : seq;
    let idx = hay.indexOf(e.s, 0);
    const seen = new Set();
    while (idx !== -1 && idx < len) {
      if (!seen.has(idx)) {
        seen.add(idx);
        sites.push({ enz: e, pos: idx + 1, cut: ((idx + e.cut) % len) + 1 });
      }
      idx = hay.indexOf(e.s, idx + 1);
    }
  }
  sites.sort((a, b) => a.pos - b.pos);
  return sites;
}

/* ---------------- 状态 ---------------- */
const state = {
  records: [],
  active: 0,
  tab: 'overview',
  strandFilter: 'both',
  topology: 'linear',
  orfs: [],
  sites: [],
  selectedOrfId: null,
  alignAuto: true,
  alignResult: null,
  mapHits: [],
};

function activeRec() { return state.records[state.active] || null; }

/* ---------------- 序列基础运算 ---------------- */
const COMP = { A: 'T', T: 'A', G: 'C', C: 'G', N: 'N', R: 'Y', Y: 'R', K: 'M', M: 'K', S: 'S', W: 'W', B: 'V', V: 'B', D: 'H', H: 'D' };
function revComp(s) {
  const out = new Array(s.length);
  for (let i = 0; i < s.length; i++) out[s.length - 1 - i] = COMP[s[i]] || 'N';
  return out.join('');
}

function parseFasta(text) {
  const recs = [];
  let cur = null;
  for (const raw of text.split(/\r?\n/)) {
    const line = raw.trim();
    if (!line) continue;
    if (line.startsWith('>')) {
      const body = line.slice(1).trim() || 'sequence_' + (recs.length + 1);
      const m = body.match(/^(\S+)\s*(.*)$/);
      cur = { name: m ? m[1] : body, desc: m ? m[2] : '', seq: '' };
      recs.push(cur);
    } else {
      if (!cur) { cur = { name: 'sequence_1', desc: '', seq: '' }; recs.push(cur); }
      cur.seq += line.replace(/[\s\d]/g, '');
    }
  }
  for (const r of recs) {
    r.seq = r.seq.toUpperCase().replace(/U/g, 'T').replace(/[^ACGTNRYKMSWBDHV]/g, '');
    if (!r.topology) r.topology = 'linear';
  }
  return recs.filter((r) => r.seq.length > 0);
}

function seqStats(seq) {
  const c = { A: 0, C: 0, G: 0, T: 0, N: 0 };
  let amb = 0;
  for (let i = 0; i < seq.length; i++) {
    const ch = seq[i];
    if (c[ch] !== undefined) c[ch]++;
    else amb++;
  }
  const acgt = c.A + c.C + c.G + c.T;
  return { len: seq.length, c, amb, acgt, gc: acgt ? (c.G + c.C) / acgt : 0 };
}

function gcOf(s) {
  let g = 0, c = 0, n = 0;
  for (let i = 0; i < s.length; i++) {
    const ch = s[i];
    if (ch === 'G') { g++; n++; } else if (ch === 'C') { c++; n++; } else if (ch === 'A' || ch === 'T') n++;
  }
  return n ? (g + c) / n : 0;
}

function gcSeries(seq, win) {
  const n = seq.length;
  win = Math.min(win, n);
  const step = Math.max(1, Math.floor(win / 4));
  const pts = [];
  for (let s = 0; s + win <= n; s += step) {
    let g = 0, c = 0;
    for (let i = s; i < s + win; i++) {
      const ch = seq[i];
      if (ch === 'G') g++; else if (ch === 'C') c++;
    }
    const gc = (g + c) / win;
    const skew = g + c ? (g - c) / (g + c) : 0;
    pts.push({ pos: s + win / 2, gc, skew });
  }
  if (pts.length > 1500) {
    const stride = Math.ceil(pts.length / 1500);
    const thin = [];
    for (let i = 0; i < pts.length; i += stride) thin.push(pts[i]);
    if (thin[thin.length - 1] !== pts[pts.length - 1]) thin.push(pts[pts.length - 1]);
    return thin;
  }
  return pts;
}

/* ---------------- ORF 预测 ---------------- */
function findORFs(seq, minAA) {
  const out = [];
  const len = seq.length;
  const minCodons = Math.max(1, minAA);
  const sides = [
    { strand: '+', s: seq },
    { strand: '-', s: revComp(seq) },
  ];
  for (const side of sides) {
    const s = side.s;
    for (let frame = 0; frame < 3; frame++) {
      const K = Math.floor((s.length - frame) / 3);
      if (K <= 0) continue;
      const aa = new Array(K);
      for (let k = 0; k < K; k++) {
        const p = frame + k * 3;
        aa[k] = CODON_AA[s.slice(p, p + 3)] || 'X';
      }
      const nextStop = new Int32Array(K + 1);
      nextStop[K] = K;
      for (let k = K - 1; k >= 0; k--) nextStop[k] = aa[k] === '*' ? k : nextStop[k + 1];

      for (let k = 0; k < K; k++) {
        if (aa[k] !== 'M') continue;
        const stop = nextStop[k + 1];
        if (stop >= K) continue;                 // 无终止子
        const nCod = stop - k;                   // 不含终止密码子
        if (nCod < minCodons) continue;
        const pStart = frame + k * 3;            // 该链 0-based
        const pEnd = frame + stop * 3 + 3;       // 含终止密码子，exclusive
        const nt = s.slice(pStart, pEnd);
        let start1, end1;
        if (side.strand === '+') { start1 = pStart + 1; end1 = pEnd; }
        else { start1 = len - pEnd + 1; end1 = len - pStart; }
        const protein = nt.slice(0, nCod * 3).match(/.{3}/g).map((c) => CODON_AA[c]).join('');
        let g3 = 0, n3 = 0;
        for (let i = 2; i < nCod * 3; i += 3) {
          const ch = nt[i];
          if (ch === 'G' || ch === 'C') g3++;
          n3++;
        }
        out.push({
          id: out.length,
          strand: side.strand,
          frame: frame + 1,
          start1, end1,
          bp: pEnd - pStart,
          aa: nCod,
          nt,
          protein,
          gc3: n3 ? g3 / n3 : 0,
        });
      }
    }
  }
  out.sort((a, b) => b.bp - a.bp);
  out.forEach((o, i) => { o.id = i; });
  return out;
}

/* ---------------- 绘制：概览 ---------------- */
function drawComposition(ctx, w, h, st) {
  const cx = w * 0.42, cy = h * 0.5;
  const R = Math.min(w * 0.32, h * 0.36);
  const inner = R * 0.62;
  const vals = [
    { k: 'A', v: st.c.A, color: '#16a34a' },
    { k: 'T', v: st.c.T, color: '#dc2626' },
    { k: 'G', v: st.c.G, color: '#f59e0b' },
    { k: 'C', v: st.c.C, color: '#2563eb' },
  ];
  const total = vals.reduce((s, x) => s + x.v, 0) || 1;
  let a0 = -Math.PI / 2;
  const gap = 0.02;
  for (const seg of vals) {
    const ang = (seg.v / total) * Math.PI * 2;
    ctx.beginPath();
    ctx.arc(cx, cy, R, a0 + gap / 2, a0 + ang - gap / 2);
    ctx.arc(cx, cy, inner, a0 + ang - gap / 2, a0 + gap / 2, true);
    ctx.closePath();
    ctx.fillStyle = seg.color;
    ctx.globalAlpha = 0.92;
    ctx.fill();
    ctx.globalAlpha = 1;
    a0 += ang;
  }
  ctx.textAlign = 'center';
  ctx.fillStyle = '#8b9aaa';
  ctx.font = '11px "Segoe UI", system-ui, sans-serif';
  ctx.fillText('GC 含量', cx, cy - 10);
  ctx.fillStyle = '#0e1c26';
  ctx.font = '600 22px "Cascadia Mono", Consolas, monospace';
  ctx.fillText((st.gc * 100).toFixed(1) + '%', cx, cy + 14);

  const lx = cx + R + 26;
  let ly = cy - 46;
  ctx.textAlign = 'left';
  for (const seg of vals) {
    ctx.fillStyle = seg.color;
    if (ctx.roundRect) { ctx.beginPath(); ctx.roundRect(lx, ly - 7, 9, 9, 2); ctx.fill(); }
    else ctx.fillRect(lx, ly - 7, 9, 9);
    ctx.fillStyle = '#516476';
    ctx.font = '11px "Segoe UI", system-ui, sans-serif';
    ctx.fillText(seg.k, lx + 15, ly + 1);
    ctx.fillStyle = '#0e1c26';
    ctx.font = '600 11px "Cascadia Mono", Consolas, monospace';
    ctx.fillText(fmt(seg.v), lx + 32, ly + 1);
    ctx.fillStyle = '#8b9aaa';
    ctx.font = '10.5px "Cascadia Mono", Consolas, monospace';
    ctx.fillText(((seg.v / total) * 100).toFixed(1) + '%', lx + 32, ly + 14);
    ly += 34;
  }
  if (st.amb) {
    ctx.fillStyle = '#8b9aaa';
    ctx.font = '10px "Segoe UI", system-ui, sans-serif';
    ctx.fillText(`含 ${st.amb} 个简并碱基`, lx, ly + 6);
  }
}

function drawGCChart(ctx, w, h, series, meanGC) {
  if (!series.length) return;
  const padL = 34, padR = 10, padT = 10, padB = 18;
  const innerW = w - padL - padR;
  const topH = Math.round((h - padT - padB) * 0.64);
  const gapY = 16;
  const botTop = padT + topH + gapY;
  const botH = h - padB - botTop;

  const n = series.length;
  const X = (i) => padL + (n <= 1 ? innerW / 2 : (i / (n - 1)) * innerW);

  ctx.strokeStyle = '#eef2f5';
  ctx.lineWidth = 1;
  for (let k = 0; k <= 4; k++) {
    const y = padT + (topH / 4) * k;
    ctx.beginPath(); ctx.moveTo(padL, y); ctx.lineTo(w - padR, y); ctx.stroke();
    ctx.fillStyle = '#8b9aaa';
    ctx.font = '9.5px "Cascadia Mono", Consolas, monospace';
    ctx.textAlign = 'right';
    ctx.fillText(String(100 - k * 25), padL - 6, y + 3);
  }

  const yOf = (v) => padT + topH - clamp(v, 0, 1) * topH;
  ctx.save();
  ctx.setLineDash([4, 4]);
  ctx.strokeStyle = '#94a3b8';
  ctx.beginPath(); ctx.moveTo(padL, yOf(meanGC)); ctx.lineTo(w - padR, yOf(meanGC)); ctx.stroke();
  ctx.restore();

  const grad = ctx.createLinearGradient(0, padT, 0, padT + topH);
  grad.addColorStop(0, 'rgba(15,118,110,.22)');
  grad.addColorStop(1, 'rgba(15,118,110,.02)');
  ctx.beginPath();
  ctx.moveTo(X(0), yOf(series[0].gc));
  for (let i = 1; i < n; i++) ctx.lineTo(X(i), yOf(series[i].gc));
  ctx.lineTo(X(n - 1), padT + topH);
  ctx.lineTo(X(0), padT + topH);
  ctx.closePath();
  ctx.fillStyle = grad;
  ctx.fill();

  ctx.beginPath();
  ctx.moveTo(X(0), yOf(series[0].gc));
  for (let i = 1; i < n; i++) ctx.lineTo(X(i), yOf(series[i].gc));
  ctx.strokeStyle = '#0f766e';
  ctx.lineWidth = 1.6;
  ctx.lineJoin = 'round';
  ctx.stroke();

  ctx.fillStyle = '#8b9aaa';
  ctx.font = '9.5px "Cascadia Mono", Consolas, monospace';
  ctx.textAlign = 'left';
  ctx.fillText('GC%', padL + 2, padT + 9);

  const zeroY = botTop + botH / 2;
  ctx.strokeStyle = '#eef2f5';
  ctx.beginPath(); ctx.moveTo(padL, zeroY); ctx.lineTo(w - padR, zeroY); ctx.stroke();
  ctx.fillStyle = '#b6c2cd';
  ctx.font = '9.5px "Cascadia Mono", Consolas, monospace';
  ctx.textAlign = 'right';
  ctx.fillText('0', padL - 6, zeroY + 3);
  ctx.textAlign = 'left';
  ctx.fillText('GC skew', padL + 2, botTop + 9);

  const ySk = (v) => zeroY - clamp(v, -1, 1) * (botH / 2 - 3);
  ctx.beginPath();
  for (let i = 0; i < n; i++) {
    const x = X(i), y = ySk(series[i].skew);
    if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
  }
  ctx.strokeStyle = '#f59e0b';
  ctx.lineWidth = 1.3;
  ctx.stroke();
  ctx.lineTo(X(n - 1), zeroY);
  ctx.lineTo(X(0), zeroY);
  ctx.closePath();
  ctx.fillStyle = 'rgba(245,158,11,.10)';
  ctx.fill();

  ctx.fillStyle = '#8b9aaa';
  ctx.font = '9.5px "Cascadia Mono", Consolas, monospace';
  ctx.textAlign = 'center';
  const total = series[n - 1].pos * 2;
  const step = niceStep(total, 6);
  for (let p = 0; p <= total; p += step) {
    const x = padL + (p / total) * innerW;
    if (x > w - padR) break;
    ctx.fillText(p >= 1000 ? (p / 1000).toFixed(p % 1000 ? 1 : 0) + 'k' : String(p), x, h - 4);
  }
}

function drawStrip(ctx, w, h, seq) {
  const len = seq.length;
  const colors = { A: '#16a34a', C: '#2563eb', G: '#f59e0b', T: '#dc2626', N: '#cbd5e1' };
  const barTop = 16, barH = h - 16;
  for (let x = 0; x < w; x++) {
    const i = Math.min(len - 1, Math.floor((x / w) * len));
    ctx.fillStyle = colors[seq[i]] || '#e2e8f0';
    ctx.globalAlpha = 0.92;
    ctx.fillRect(x, barTop, 1, barH);
  }
  ctx.globalAlpha = 1;
  ctx.fillStyle = '#8b9aaa';
  ctx.font = '9.5px "Cascadia Mono", Consolas, monospace';
  const step = niceStep(len, Math.max(4, Math.floor(w / 90)));
  for (let p = 0; p <= len; p += step) {
    const x = (p / len) * w;
    if (x > w - 2) break;
    ctx.fillRect(x, barTop - 6, 1, 4);
    ctx.textAlign = x < 20 ? 'left' : x > w - 40 ? 'right' : 'center';
    ctx.fillText(p >= 1000 ? (p / 1000).toFixed(p % 1000 ? 1 : 0) + 'k' : String(p), x, barTop - 9);
  }
}

/* ---------------- ORF 地图 ---------------- */
function renderOrfMap() {
  const rec = activeRec();
  const box = $('#orfMap');
  if (!rec) { box.innerHTML = ''; return; }
  const len = rec.seq.length;
  const orfs = state.orfs;
  const lanes = [];
  if (state.strandFilter !== 'minus') lanes.push('+1', '+2', '+3');
  lanes.push('ruler');
  if (state.strandFilter !== 'plus') lanes.push('-1', '-2', '-3');

  const frag = document.createDocumentFragment();
  for (const lane of lanes) {
    const row = document.createElement('div');
    row.className = 'lane' + (lane === 'ruler' ? ' ruler' : '');
    const lab = document.createElement('span');
    lab.className = 'lane-label';
    lab.textContent = lane === 'ruler' ? '' : (lane[0] === '+' ? '+' + lane[1] : '−' + lane[1]);
    const track = document.createElement('div');
    track.className = 'lane-track';

    if (lane === 'ruler') {
      const step = niceStep(len, 10);
      for (let p = 0; p <= len; p += step) {
        const t = document.createElement('div');
        t.className = 'ruler-tick';
        t.style.left = (p / len) * 100 + '%';
        const b = document.createElement('b');
        b.textContent = p >= 1000 ? (p / 1000).toFixed(p % 1000 ? 1 : 0) + 'kb' : p + 'bp';
        t.appendChild(b);
        track.appendChild(t);
      }
    } else {
      const strand = lane[0];
      const frame = Number(lane[1]);
      for (const o of orfs) {
        if (o.strand !== strand || o.frame !== frame) continue;
        const bar = document.createElement('div');
        bar.className = 'orf-bar ' + (o.strand === '+' ? 'plus' : 'minus');
        bar.style.left = ((o.start1 - 1) / len) * 100 + '%';
        bar.style.width = Math.max(0.35, (o.bp / len) * 100) + '%';
        bar.dataset.id = o.id;
        if (state.selectedOrfId === o.id) bar.classList.add('sel');
        bar.addEventListener('mousemove', (ev) => {
          showTip(`<b>ORF #${o.id + 1}</b> · ${o.strand}${o.frame} 框<br>
            位置 <span class="tt-mono">${fmt(o.start1)}–${fmt(o.end1)}</span><br>
            长度 <span class="tt-mono">${fmt(o.bp)} bp / ${o.aa} aa</span>`, ev.clientX, ev.clientY);
        });
        bar.addEventListener('mouseleave', hideTip);
        bar.addEventListener('click', () => selectOrf(o.id));
        track.appendChild(bar);
      }
    }
    row.appendChild(lab);
    row.appendChild(track);
    frag.appendChild(row);
  }
  box.innerHTML = '';
  box.appendChild(frag);

  const plus = orfs.filter((o) => o.strand === '+').length;
  const minus = orfs.length - plus;
  $('#orfSummary').textContent = `共 ${orfs.length} 个 ORF · 正链 ${plus} · 负链 ${minus} · 最短 ${$('#minAA').value} aa`;
  $('#orfListCount').textContent = `(${Math.min(orfs.length, 80)}/${orfs.length})`;
  renderOrfTable();
}

function renderOrfTable() {
  const tbody = $('#orfTable tbody');
  const rows = state.orfs.slice(0, 80);
  tbody.innerHTML = rows.map((o) => `
    <tr data-id="${o.id}" class="${state.selectedOrfId === o.id ? 'sel' : ''}">
      <td class="mono">${o.id + 1}</td>
      <td><span class="pill ${o.strand === '+' ? 'plus' : 'minus'}">${o.strand}${o.frame}</span></td>
      <td class="mono">${fmt(o.start1)}</td>
      <td class="mono">${fmt(o.end1)}</td>
      <td class="mono">${fmt(o.bp)} bp</td>
      <td class="mono">${o.aa} aa</td>
    </tr>`).join('');
  $$('#orfTable tbody tr').forEach((tr) => {
    tr.addEventListener('click', () => selectOrf(Number(tr.dataset.id)));
  });
}

function selectOrf(id) {
  state.selectedOrfId = id;
  const o = state.orfs.find((x) => x.id === id);
  if (!o) return;
  $$('.orf-bar').forEach((b) => b.classList.toggle('sel', Number(b.dataset.id) === id));
  $$('#orfTable tbody tr').forEach((tr) => tr.classList.toggle('sel', Number(tr.dataset.id) === id));

  const rec = activeRec();
  const nt = o.nt;
  const nCod = Math.floor(nt.length / 3);
  const shown = Math.min(nCod, 260);
  const chips = [];
  for (let i = 0; i < shown; i++) {
    const codon = nt.slice(i * 3, i * 3 + 3);
    const aa = CODON_AA[codon] || 'X';
    const cls = i === 0 ? 'start' : aa === '*' ? 'stop' : '';
    chips.push(`<span class="cod ${cls}"><span class="nt">${codon}</span><span class="aa">${aa === '*' ? '■' : aa}</span></span>`);
  }
  $('#orfMeta').textContent = `${o.strand}${o.frame} 框 · ${fmt(o.start1)}–${fmt(o.end1)} · ${o.aa} aa`;
  $('#orfDetail').innerHTML = `
    <div class="orf-head">
      <span class="chip">长度 <b>${fmt(o.bp)} bp</b></span>
      <span class="chip">氨基酸 <b>${o.aa} aa</b></span>
      <span class="chip">GC3s <b>${pctS(o.gc3)}</b></span>
      <span class="chip">起始 <b>ATG</b></span>
      <span class="chip">终止 <b>${o.nt.slice(-3)}</b></span>
    </div>
    <div class="orf-seq">${chips.join('')}${nCod > shown ? `<span class="hint" style="align-self:center;margin-left:8px">… 仅显示前 ${shown} 个密码子</span>` : ''}</div>
    <div class="row gap" style="margin-top:10px">
      <button class="btn ghost" id="btnCopyProt">复制蛋白序列</button>
      <button class="btn ghost" id="btnCopyNt">复制核苷酸序列</button>
      <span class="hint">${rec ? rec.name : ''}</span>
    </div>`;
  $('#btnCopyProt').addEventListener('click', () => copyText(o.protein));
  $('#btnCopyNt').addEventListener('click', () => copyText(o.nt));
}

function copyText(t) {
  navigator.clipboard?.writeText(t).then(
    () => showToast('已复制到剪贴板'),
    () => showToast('复制失败，请手动选择')
  );
}

let toastTimer = null;
function showToast(msg) {
  let box = $('#toast');
  if (!box) {
    box = document.createElement('div');
    box.id = 'toast';
    box.className = 'toast';
    document.body.appendChild(box);
  }
  box.textContent = msg;
  box.classList.add('show');
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => box.classList.remove('show'), 1600);
}

/* ---------------- 酶切图谱 ---------------- */
const GC_STOPS = [
  [0, [191, 219, 254]],
  [0.35, [96, 165, 250]],
  [0.5, [45, 212, 191]],
  [0.65, [5, 150, 105]],
  [1, [6, 78, 59]],
];
const CODON_STOPS = [
  [0, [238, 242, 244]],
  [0.35, [153, 230, 218]],
  [0.7, [15, 118, 110]],
  [1, [11, 79, 74]],
];
function rampColor(stops, v) {
  v = clamp(v, 0, 1);
  let a = stops[0], b = stops[stops.length - 1];
  for (let i = 0; i < stops.length - 1; i++) {
    if (v >= stops[i][0] && v <= stops[i + 1][0]) { a = stops[i]; b = stops[i + 1]; break; }
  }
  const t = (v - a[0]) / (b[0] - a[0] || 1);
  const c = a[1].map((x, i) => Math.round(x + (b[1][i] - x) * t));
  return `rgb(${c[0]},${c[1]},${c[2]})`;
}
const gcColor = (v) => rampColor(GC_STOPS, v);

function renderEnzyme() {
  const rec = activeRec();
  if (!rec) return;
  const canvas = $('#mapCanvas');
  const { ctx, w, h } = canvasCtx(canvas);
  const len = rec.seq.length;
  const sites = state.sites;
  state.mapHits = [];

  if (state.topology === 'circular') {
    const cx = w / 2, cy = h / 2;
    const R = Math.max(60, Math.min(w, h) / 2 - 62);
    const innerR = R - 42;

    ctx.beginPath();
    ctx.arc(cx, cy, R, 0, Math.PI * 2);
    ctx.strokeStyle = '#eef3f6';
    ctx.lineWidth = 18;
    ctx.stroke();
    ctx.beginPath();
    ctx.arc(cx, cy, R, 0, Math.PI * 2);
    ctx.strokeStyle = '#dfe8ee';
    ctx.lineWidth = 1;
    ctx.stroke();

    const slices = 240;
    for (let i = 0; i < slices; i++) {
      const s = Math.floor((i / slices) * len);
      const e = Math.max(s + 1, Math.floor(((i + 1) / slices) * len));
      const gc = gcOf(rec.seq.slice(s, e));
      const a0 = (i / slices) * Math.PI * 2 - Math.PI / 2;
      const a1 = ((i + 1) / slices) * Math.PI * 2 - Math.PI / 2;
      ctx.beginPath();
      ctx.arc(cx, cy, innerR, a0, a1 + 0.004);
      ctx.strokeStyle = gcColor(gc);
      ctx.lineWidth = 13;
      ctx.stroke();
    }
    ctx.beginPath();
    ctx.arc(cx, cy, innerR, 0, Math.PI * 2);
    ctx.strokeStyle = 'rgba(255,255,255,.9)';
    ctx.lineWidth = 1;
    ctx.stroke();

    const step = niceStep(len, 12);
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    for (let p = 0; p < len; p += step) {
      const ang = (p / len) * Math.PI * 2 - Math.PI / 2;
      const c = Math.cos(ang), s = Math.sin(ang);
      ctx.beginPath();
      ctx.moveTo(cx + c * (R + 9), cy + s * (R + 9));
      ctx.lineTo(cx + c * (R + 15), cy + s * (R + 15));
      ctx.strokeStyle = '#c8d4dd';
      ctx.lineWidth = 1;
      ctx.stroke();
      ctx.fillStyle = '#93a3b2';
      ctx.font = '9px "Cascadia Mono", Consolas, monospace';
      ctx.fillText(p >= 1000 ? (p / 1000).toFixed(p % 1000 ? 1 : 0) + 'k' : String(p), cx + c * (R + 25), cy + s * (R + 25));
    }
    ctx.textBaseline = 'alphabetic';

    const counts = {};
    for (const s of sites) counts[s.enz.n] = (counts[s.enz.n] || 0) + 1;
    const labelMode = $('#labelMode').value;
    const labelSet = new Set();
    if (labelMode === 'all') sites.forEach((s) => labelSet.add(s));
    else if (labelMode === 'unique') sites.forEach((s) => { if (counts[s.enz.n] === 1) labelSet.add(s); });

    for (const s of sites) {
      const ang = ((s.pos - 1) / len) * Math.PI * 2 - Math.PI / 2;
      const c = Math.cos(ang), sn = Math.sin(ang);
      ctx.beginPath();
      ctx.moveTo(cx + c * (R - 10), cy + sn * (R - 10));
      ctx.lineTo(cx + c * (R + 12), cy + sn * (R + 12));
      ctx.strokeStyle = s.enz.four ? '#94a3b8' : '#0f766e';
      ctx.lineWidth = s.enz.four ? 1 : 1.8;
      ctx.stroke();
      state.mapHits.push({ x: cx + c * (R + 12), y: cy + sn * (R + 12), r: 9, site: s });
    }

    /* 标签排布：密集的位点自动聚合为区间标注，稀疏位点分层外挂 */
    const items = sites
      .filter((s) => labelSet.has(s))
      .map((s) => ({ s, ang: ((s.pos - 1) / len) * Math.PI * 2 - Math.PI / 2 }))
      .sort((a, b) => a.ang - b.ang);

    const groups = [];
    for (const it of items) {
      const g = groups[groups.length - 1];
      if (g && it.ang - g[g.length - 1].ang < 0.13) g.push(it);
      else groups.push([it]);
    }

    const drawRadialLabel = (it, radius) => {
      const { ang } = it;
      const c = Math.cos(ang), sn = Math.sin(ang);
      ctx.save();
      ctx.translate(cx + c * radius, cy + sn * radius);
      const flip = c < 0;
      ctx.rotate(ang + (flip ? Math.PI : 0));
      ctx.textAlign = flip ? 'right' : 'left';
      ctx.textBaseline = 'middle';
      ctx.fillStyle = '#33475b';
      ctx.font = '9.5px "Cascadia Mono", Consolas, monospace';
      ctx.fillText(it.s.enz.n, 0, 0);
      ctx.restore();
    };

    for (const g of groups) {
      if (g.length >= 4) {
        const a0 = g[0].ang, a1 = g[g.length - 1].ang;
        const rB = R + 22;
        ctx.beginPath();
        ctx.arc(cx, cy, rB, a0, a1);
        ctx.strokeStyle = '#9fb0bf';
        ctx.lineWidth = 1.1;
        ctx.stroke();
        for (const a of [a0, a1]) {
          ctx.beginPath();
          ctx.moveTo(cx + Math.cos(a) * (R + 10), cy + Math.sin(a) * (R + 10));
          ctx.lineTo(cx + Math.cos(a) * (rB + 2), cy + Math.sin(a) * (rB + 2));
          ctx.strokeStyle = '#9fb0bf';
          ctx.lineWidth = 1.1;
          ctx.stroke();
        }
        const mid = (a0 + a1) / 2;
        const lx = cx + Math.cos(mid) * (rB + 20);
        const ly = cy + Math.sin(mid) * (rB + 20);
        const text = `${g.length} 个位点`;
        ctx.font = '600 10px "Segoe UI", system-ui, sans-serif';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.lineWidth = 3;
        ctx.strokeStyle = 'rgba(255,255,255,.92)';
        ctx.strokeText(text, lx, ly);
        ctx.fillStyle = '#0f766e';
        ctx.fillText(text, lx, ly);
      } else {
        const lastAt = [];
        for (const it of g) {
          let lv = 0;
          while (lv < 3 && lastAt[lv] !== undefined && it.ang - lastAt[lv] < 0.13) lv++;
          if (lv >= 3) continue;
          lastAt[lv] = it.ang;
          const radius = R + 30 + lv * 15;
          if (lv > 0) {
            ctx.beginPath();
            ctx.moveTo(cx + Math.cos(it.ang) * (R + 12), cy + Math.sin(it.ang) * (R + 12));
            ctx.lineTo(cx + Math.cos(it.ang) * (radius - 2), cy + Math.sin(it.ang) * (radius - 2));
            ctx.strokeStyle = '#cbd5e1';
            ctx.lineWidth = 0.9;
            ctx.stroke();
          }
          drawRadialLabel(it, radius);
        }
      }
    }

    ctx.textAlign = 'center';
    ctx.textBaseline = 'alphabetic';
    ctx.fillStyle = '#0e1c26';
    ctx.font = '600 15px "Segoe UI", system-ui, sans-serif';
    ctx.fillText(rec.name, cx, cy - 14);
    ctx.fillStyle = '#516476';
    ctx.font = '12px "Cascadia Mono", Consolas, monospace';
    ctx.fillText(`${fmt(len)} bp · 环状`, cx, cy + 6);
    ctx.fillStyle = '#8b9aaa';
    ctx.font = '10.5px "Segoe UI", system-ui, sans-serif';
    const nEnz = Object.keys(counts).length;
    ctx.fillText(`${nEnz} 种酶 / ${sites.length} 个位点`, cx, cy + 24);
    const uniq = Object.keys(counts).filter((k) => counts[k] === 1).length;
    if (uniq) ctx.fillText(`单切点酶 ${uniq} 种`, cx, cy + 40);
  } else {
    drawLinearMap(ctx, w, h, rec, sites);
  }
}

function drawLinearMap(ctx, w, h, rec, sites) {
  const len = rec.seq.length;
  const x0 = 52, x1 = w - 46;
  const y = Math.round(h * 0.46);
  const X = (p) => x0 + (p / len) * (x1 - x0);

  ctx.fillStyle = '#e8eef3';
  ctx.fillRect(x0, y - 5, x1 - x0, 10);
  ctx.fillStyle = '#d5e0e8';
  ctx.fillRect(x0, y - 5, x1 - x0, 1);

  const step = niceStep(len, 12);
  ctx.textAlign = 'center';
  ctx.textBaseline = 'alphabetic';
  for (let p = 0; p <= len; p += step) {
    const x = X(p);
    ctx.fillStyle = '#c8d4dd';
    ctx.fillRect(x, y + 6, 1, 7);
    ctx.fillStyle = '#93a3b2';
    ctx.font = '9px "Cascadia Mono", Consolas, monospace';
    ctx.fillText(p >= 1000 ? (p / 1000).toFixed(p % 1000 ? 1 : 0) + 'k' : String(p), x, y + 25);
  }

  const counts = {};
  for (const s of sites) counts[s.enz.n] = (counts[s.enz.n] || 0) + 1;
  const labelMode = $('#labelMode').value;
  let row = 0;
  for (const s of sites) {
    const x = X(s.cut - 1);
    const up = row % 2 === 0;
    ctx.beginPath();
    ctx.moveTo(x, y - 5);
    ctx.lineTo(x, up ? y - 22 : y + 22);
    ctx.strokeStyle = s.enz.four ? '#cbd5e1' : '#0f766e';
    ctx.lineWidth = s.enz.four ? 1 : 1.5;
    ctx.stroke();
    state.mapHits.push({ x, y: up ? y - 22 : y + 22, r: 8, site: s });
    if (labelMode === 'all' || (labelMode === 'unique' && counts[s.enz.n] === 1)) {
      ctx.save();
      ctx.translate(x, up ? y - 26 : y + 26);
      ctx.rotate(-Math.PI / 3);
      ctx.textAlign = up ? 'left' : 'right';
      ctx.fillStyle = '#33475b';
      ctx.font = '10px "Cascadia Mono", Consolas, monospace';
      ctx.fillText(s.enz.n, 0, 0);
      ctx.restore();
    }
    row++;
  }

  ctx.textAlign = 'left';
  ctx.fillStyle = '#0e1c26';
  ctx.font = '600 13px "Segoe UI", system-ui, sans-serif';
  ctx.fillText(rec.name, x0, 22);
  ctx.fillStyle = '#8b9aaa';
  ctx.font = '10.5px "Cascadia Mono", Consolas, monospace';
  ctx.fillText(`${fmt(len)} bp · 线性`, x0 + 4, 38);
}

function renderSiteTable() {
  const tbody = $('#siteTable tbody');
  const sites = state.sites;
  $('#siteCount').textContent = sites.length + ' 个';
  tbody.innerHTML = sites.slice(0, 400).map((s) => {
    const end = enzymeEnd(s.enz);
    return `<tr>
      <td><b>${s.enz.n}</b></td>
      <td class="mono">${s.enz.s}</td>
      <td class="mono">${fmt(s.pos)}</td>
      <td class="mono">${fmt(s.cut)}</td>
      <td><span class="pill ${end.cls}">${end.t}</span></td>
    </tr>`;
  }).join('');
}

/* ---------------- 密码子 ---------------- */
function codonStats(orf) {
  const nt = orf.nt;
  const counts = Object.create(null);
  let total = 0;
  for (let i = 0; i + 3 <= nt.length; i += 3) {
    const c = nt.slice(i, i + 3);
    const aa = CODON_AA[c];
    if (!aa || aa === '*') continue;
    counts[c] = (counts[c] || 0) + 1;
    total++;
  }
  const aaCount = Object.create(null);
  for (const c in counts) aaCount[CODON_AA[c]] = (aaCount[CODON_AA[c]] || 0) + counts[c];

  const rscu = Object.create(null);
  for (const aa in SYN) {
    if (aa === '*') continue;
    const syn = SYN[aa];
    const tot = aaCount[aa] || 0;
    for (const c of syn) rscu[c] = tot ? ((counts[c] || 0) * syn.length) / tot : 0;
  }

  let gc3 = 0, n3 = 0;
  for (const c in counts) {
    const syn = SYN[CODON_AA[c]];
    if (syn.length === 1) continue;
    const b = c[2];
    if (b === 'G' || b === 'C') gc3 += counts[c];
    n3 += counts[c];
  }
  gc3 = n3 ? gc3 / n3 : 0;

  const F = { 2: [], 3: [], 4: [], 6: [] };
  for (const aa in SYN) {
    if (aa === '*') continue;
    const syn = SYN[aa];
    const k = syn.length;
    if (!F[k]) continue;
    const n = aaCount[aa] || 0;
    if (n < 2) continue;
    let sumP2 = 0;
    for (const c of syn) { const p = (counts[c] || 0) / n; sumP2 += p * p; }
    F[k].push((n * sumP2 - 1) / (n - 1));
  }
  const avg = (arr) => (arr.length ? arr.reduce((a, b) => a + b, 0) / arr.length : null);
  const f2 = avg(F[2]), f3 = avg(F[3]), f4 = avg(F[4]), f6 = avg(F[6]);
  let enc = 2;
  if (f2) enc += 9 / f2;
  if (f3) enc += 1 / f3;
  if (f4) enc += 5 / f4;
  if (f6) enc += 3 / f6;
  enc = Math.min(enc, 61);

  return { counts, total, aaCount, rscu, gc3, enc };
}

function renderCodon() {
  const sel = $('#orfPicker');
  const orf = state.orfs[Number(sel.value)] || state.orfs[0];
  if (!orf) return;
  const cs = codonStats(orf);
  $('#cdTotal').textContent = fmt(cs.total);
  $('#cdGC3').textContent = pctS(cs.gc3);
  $('#cdENC').textContent = cs.enc.toFixed(1);

  const parts = [];
  for (let r = 0; r < 16; r++) {
    const p1 = BASE_ORDER[Math.floor(r / 4)];
    const p2 = BASE_ORDER[r % 4];
    parts.push(`<div class="codon-prefix">${p1}${p2}</div>`);
    for (let k = 0; k < 4; k++) {
      const codon = p1 + p2 + BASE_ORDER[k];
      const aa = CODON_AA[codon];
      if (aa === '*') {
        parts.push(`<div class="codon-cell stop"><span class="cod-name">${codon}</span><span class="cod-val">Stop</span></div>`);
        continue;
      }
      const n = cs.counts[codon] || 0;
      const syn = SYN[aa].length;
      const v = clamp(cs.rscu[codon] / syn, 0, 1);
      const bg = rampColor(CODON_STOPS, v);
      const rgb = bg.match(/\d+/g).map(Number);
      const lum = (0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]) / 255;
      const dark = lum < 0.62;
      parts.push(`<div class="codon-cell" style="background:${bg};color:${dark ? '#f8fafc' : '#0e1c26'}"
        title="${aa} (${AA_NAMES[aa]}) · ${n} 次 · RSCU ${cs.rscu[codon].toFixed(2)}">
        <span class="cod-name">${codon}</span>
        <span class="cod-val">${n} · ${cs.rscu[codon].toFixed(2)}</span></div>`);
    }
  }
  $('#codonGrid').innerHTML = parts.join('');

  const bars = AA20.map((aa) => ({ aa, n: cs.aaCount[aa] || 0 })).sort((a, b) => b.n - a.n);
  const max = Math.max(1, ...bars.map((b) => b.n));
  $('#aaBars').innerHTML = bars.map((b) => `
    <div class="aa-row">
      <span class="aa-name">${b.aa}</span>
      <span class="aa-bar"><i style="width:${(b.n / max) * 100}%"></i></span>
      <span class="aa-val">${fmt(b.n)} · ${pctS(b.n / Math.max(1, cs.total))}</span>
    </div>`).join('');
}

/* ---------------- 双序列比对 ---------------- */
function alignPair(a, b, opt) {
  const n = a.length, m = b.length, W = m + 1;
  const size = (n + 1) * W;
  const H = new Int32Array(size);
  const T = new Uint8Array(size); // 1=diag 2=up 3=left
  const gap = opt.gap, match = opt.match, mism = opt.mism;
  const local = opt.mode === 'local';
  let best = { score: local ? 0 : -Infinity, i: n, j: m };

  for (let j = 0; j <= m; j++) { H[j] = local ? 0 : j * gap; T[j] = 3; }
  for (let i = 1; i <= n; i++) { H[i * W] = local ? 0 : i * gap; T[i * W] = 2; }

  for (let i = 1; i <= n; i++) {
    const ai = a[i - 1];
    const row = i * W, prev = (i - 1) * W;
    for (let j = 1; j <= m; j++) {
      const d = H[prev + j - 1] + (ai === b[j - 1] ? match : mism);
      let v = d, t = 1;
      const u = H[prev + j] + gap;
      if (u > v) { v = u; t = 2; }
      const l = H[row + j - 1] + gap;
      if (l > v) { v = l; t = 3; }
      if (local && v < 0) { v = 0; t = 0; }
      H[row + j] = v;
      T[row + j] = t;
      if (local && v > best.score) best = { score: v, i, j };
    }
  }
  if (!local) best.score = H[n * W + m];

  let i = best.i, j = best.j;
  const ca = [], cb = [];
  while (i > 0 && j > 0) {
    const t = T[i * W + j];
    if (t === 0) break;
    if (t === 1) { ca.push(a[i - 1]); cb.push(b[j - 1]); i--; j--; }
    else if (t === 2) { ca.push(a[i - 1]); cb.push('-'); i--; }
    else { ca.push('-'); cb.push(b[j - 1]); j--; }
  }
  if (!local) {
    while (i > 0) { ca.push(a[i - 1]); cb.push('-'); i--; }
    while (j > 0) { ca.push('-'); cb.push(b[j - 1]); j--; }
  }
  ca.reverse(); cb.reverse();
  return { a: ca.join(''), b: cb.join(''), score: best.score, startA: i, startB: j };
}

function drawDot(seqA, seqB) {
  const canvas = $('#dotCanvas');
  const { ctx, w, h } = canvasCtx(canvas);
  const padL = 42, padB = 32, padT = 8, padR = 8;
  const W = w - padL - padR, Hh = h - padT - padB;
  const k = 8;

  ctx.fillStyle = '#fbfdfe';
  ctx.fillRect(padL, padT, W, Hh);
  ctx.strokeStyle = '#e6edf3';
  ctx.lineWidth = 1;
  ctx.strokeRect(padL, padT, W, Hh);

  ctx.strokeStyle = '#f0f4f7';
  for (let g = 1; g < 5; g++) {
    const x = padL + (W / 5) * g, y = padT + (Hh / 5) * g;
    ctx.beginPath(); ctx.moveTo(x, padT); ctx.lineTo(x, padT + Hh); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(padL, y); ctx.lineTo(padL + W, y); ctx.stroke();
  }

  if (seqA.length >= k && seqB.length >= k) {
    const idx = new Map();
    for (let i = 0; i + k <= seqA.length; i++) {
      const key = seqA.slice(i, i + k);
      let arr = idx.get(key);
      if (!arr) { arr = []; idx.set(key, arr); }
      if (arr.length < 40) arr.push(i);
    }
    for (let j = 0; j + k <= seqB.length; j++) {
      const arr = idx.get(seqB.slice(j, j + k));
      if (!arr) continue;
      for (const ii of arr) {
        const x = padL + (ii / seqA.length) * W;
        const y = padT + Hh - (j / seqB.length) * Hh;
        ctx.fillStyle = 'rgba(15,118,110,.55)';
        ctx.fillRect(x - 1, y - 1, 2.2, 2.2);
      }
    }
  }

  ctx.fillStyle = '#8b9aaa';
  ctx.font = '9.5px "Cascadia Mono", Consolas, monospace';
  ctx.textAlign = 'center';
  for (let g = 0; g <= 4; g++) {
    ctx.fillText(fmt(Math.round((seqA.length / 4) * g)), padL + (W / 4) * g, h - 12);
  }
  ctx.save();
  ctx.translate(12, padT + Hh / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = 'center';
  ctx.fillText(fmt(seqB.length), 0, 0);
  ctx.restore();
  ctx.textAlign = 'left';
  ctx.fillText('A →', padL, h - 2);
}

function renderAlign() {
  const A = $('#alnA').value.replace(/[^A-Za-z]/g, '').toUpperCase();
  const B = $('#alnB').value.replace(/[^A-Za-z]/g, '').toUpperCase();
  const view = $('#alnView');
  const statsBox = $('#alnStats');
  if (!A || !B) {
    view.innerHTML = '<p class="empty">请先在两个输入框中填入序列。</p>';
    statsBox.innerHTML = '';
    drawDot('', '');
    return;
  }
  if (A.length * B.length > 12e6) {
    view.innerHTML = `<p class="empty">序列过长（${A.length} × ${B.length}），请截取片段后再比对。</p>`;
    return;
  }
  const mode = $('#alnMode .active').dataset.mode;
  const opt = {
    mode,
    match: Number($('#sMatch').value) || 0,
    mism: Number($('#sMis').value) || 0,
    gap: Number($('#sGap').value) || 0,
  };
  const t0 = performance.now();
  const r = alignPair(A, B, opt);
  const dt = performance.now() - t0;
  state.alignResult = r;

  let ident = 0, gaps = 0;
  for (let i = 0; i < r.a.length; i++) {
    const x = r.a[i], y = r.b[i];
    if (x === '-' || y === '-') gaps++;
    else if (x === y) ident++;
  }
  const cols = r.a.length;
  const mism = cols - ident - gaps;
  statsBox.innerHTML = `
    <span class="stat-pill">比对得分 <b>${r.score}</b></span>
    <span class="stat-pill">比对长度 <b>${fmt(cols)}</b></span>
    <span class="stat-pill">一致度 <b>${cols ? ((ident / cols) * 100).toFixed(1) : 0}%</b></span>
    <span class="stat-pill">一致 / 错配 / 空位 <b>${fmt(ident)} / ${fmt(mism)} / ${fmt(gaps)}</b></span>
    <span class="stat-pill">耗时 <b>${dt.toFixed(0)} ms</b></span>`;

  const perLine = 60;
  const maxCols = 4200;
  const shown = Math.min(cols, maxCols);
  const lines = [];
  for (let s = 0; s < shown; s += perLine) {
    const e = Math.min(shown, s + perLine);
    const segA = r.a.slice(s, e), segB = r.b.slice(s, e);
    let ra = '', rb = '', rm = '';
    for (let i = 0; i < segA.length; i++) {
      const x = segA[i], y = segB[i];
      const cls = x === '-' || y === '-' ? 'g' : x === y ? 'm' : 'x';
      ra += `<span class="${cls}">${x}</span>`;
      rb += `<span class="${cls}">${y}</span>`;
      rm += `<span class="${cls}">${cls === 'm' ? '│' : ' '}</span>`;
    }
    lines.push(
      `<span class="aln-line"><span class="aln-gutter">A ${String(r.startA + s + 1).padStart(6)}</span>${ra}</span>`,
      `<span class="aln-line"><span class="aln-gutter">  ${' '.repeat(6)}</span>${rm}</span>`,
      `<span class="aln-line"><span class="aln-gutter">B ${String(r.startB + s + 1).padStart(6)}</span>${rb}</span>`,
      `<span class="aln-line">&nbsp;</span>`
    );
  }
  if (cols > shown) lines.push(`<span class="aln-line muted">… 共 ${fmt(cols)} 列，仅显示前 ${fmt(shown)} 列</span>`);
  view.innerHTML = lines.join('');
  drawDot(A, B);
}

/* ---------------- 主流程 ---------------- */
function recStats(rec) {
  if (!rec._stats) rec._stats = seqStats(rec.seq);
  return rec._stats;
}

function getOrfs(rec) {
  const minAA = Number($('#minAA').value);
  const key = rec.seq.length + ':' + minAA;
  if (!rec._orfCache || rec._orfCache.key !== key) {
    rec._orfCache = { key, list: findORFs(rec.seq, minAA) };
  }
  return rec._orfCache.list;
}

function getSites(rec) {
  const set = $('#enzymeSet').value;
  const topo = state.topology;
  const excludeFour = $('#primeOnly').checked;
  let pool = ENZYMES;
  if (set === 'common') pool = ENZYMES.filter((e) => e.c);
  const key = topo + ':' + set + ':' + excludeFour;
  if (!rec._siteCache) rec._siteCache = new Map();
  if (rec._siteCache.has(key)) return rec._siteCache.get(key);
  let list = findSites(rec.seq, pool, topo === 'circular');
  if (set === 'unique') {
    const counts = {};
    for (const s of list) counts[s.enz.n] = (counts[s.enz.n] || 0) + 1;
    list = list.filter((s) => counts[s.enz.n] === 1);
  }
  if (excludeFour) list = list.filter((s) => !s.enz.four);
  rec._siteCache.set(key, list);
  return list;
}

function refreshAll() {
  const rec = activeRec();
  if (!rec) return;
  const st = recStats(rec);
  const tooLong = rec.seq.length > 200000;

  $('#kLen').textContent = fmt(st.len);
  $('#kGC').textContent = (st.gc * 100).toFixed(1);
  $('#statusChip').classList.add('ok');
  $('#statusChip').classList.remove('warn');
  $('#statusChip').querySelector('span').textContent =
    `${rec.name} · ${fmt(st.len)} bp · ${state.topology === 'circular' ? '环状' : '线性'}`;

  if (tooLong) {
    state.orfs = [];
    state.sites = [];
    $('#orfSummary').textContent = '序列超过 200 kb，已跳过 ORF / 酶切分析';
  } else {
    state.orfs = getOrfs(rec);
    state.sites = getSites(rec);
  }
  if (!state.orfs.some((o) => o.id === state.selectedOrfId)) {
    state.selectedOrfId = state.orfs.length ? state.orfs[0].id : null;
  }

  $('#kOrf').textContent = fmt(state.orfs.length);
  $('#kSite').textContent = fmt(state.sites.length);

  // 卡片副标题与底栏
  $('#kLenSub').textContent = state.topology === 'circular' ? '环状序列' : '线性序列';
  $('#kGCSub').textContent = `滑窗 ${$('#gcWin').value} bp`;
  let longest = 0;
  for (const o of state.orfs) if (o.aa > longest) longest = o.aa;
  $('#kOrfSub').textContent = longest ? `最长 ${fmt(longest)} aa` : '未检出';
  const counts = {};
  for (const s of state.sites) counts[s.enz.n] = (counts[s.enz.n] || 0) + 1;
  const uniqEnz = Object.keys(counts).filter((k) => counts[k] === 1).length;
  $('#kSiteSub').textContent = uniqEnz ? `单切点 ${uniqEnz} 种` : `${Object.keys(counts).length} 种酶`;
  $('#footMeta').textContent = `${rec.name} · ${fmt(st.len)} bp · ${rec.desc || '本地序列'}`;

  const picker = $('#orfPicker');
  const options = state.orfs.slice(0, 40);
  picker.innerHTML = options.length
    ? options.map((o) => `<option value="${o.id}">ORF ${o.id + 1} · ${o.strand}${o.frame} · ${o.aa} aa · ${fmt(o.start1)}–${fmt(o.end1)}</option>`).join('')
    : '<option>无 ORF</option>';
  if (state.selectedOrfId != null) picker.value = String(state.selectedOrfId);

  renderOrfMap();
  renderSiteTable();
  if (state.selectedOrfId != null) selectOrf(state.selectedOrfId);
  else $('#orfDetail').innerHTML = '<p class="empty">当前参数下没有找到 ORF，可尝试降低最短长度。</p>';

  renderActiveTab();
}

function renderActiveTab() {
  const rec = activeRec();
  if (!rec) return;
  const st = recStats(rec);
  if (state.tab === 'overview') {
    const c1 = canvasCtx($('#compCanvas'));
    drawComposition(c1.ctx, c1.w, c1.h, st);
    const win = Number($('#gcWin').value);
    if (!rec._gcs || rec._gcs.win !== win) rec._gcs = { win, data: gcSeries(rec.seq, win) };
    const c2 = canvasCtx($('#gcCanvas'));
    drawGCChart(c2.ctx, c2.w, c2.h, rec._gcs.data, st.gc);
    const c3 = canvasCtx($('#stripCanvas'));
    drawStrip(c3.ctx, c3.w, c3.h, rec.seq);
    $('#gcHint').textContent = `窗口 ${win} bp`;
    $('#stripHint').textContent = `等宽压缩显示 · ${fmt(st.len)} bp`;
  } else if (state.tab === 'enzyme') {
    renderEnzyme();
  } else if (state.tab === 'codon') {
    renderCodon();
  } else if (state.tab === 'align') {
    if ($('#alnA').value && $('#alnB').value) renderAlign();
  }
}

/* ---------------- 事件 ---------------- */
function mutate(seq, seed) {
  const rnd = mulberry32(seed || 7);
  const bases = 'ACGT';
  const out = seq.split('');
  for (let i = 0; i < out.length; i++) {
    if (rnd() < 0.035) out[i] = bases[Math.floor(rnd() * 4)];
  }
  const delAt = Math.floor(rnd() * Math.max(1, out.length - 24));
  out.splice(delAt, 6);
  const insAt = Math.floor(rnd() * Math.max(1, out.length - 12));
  const ins = Array.from({ length: 8 }, () => bases[Math.floor(rnd() * 4)]);
  out.splice(insAt, 0, ...ins);
  return out.join('');
}

function applyTopologyButtons() {
  $$('#topoSeg button').forEach((b) => b.classList.toggle('active', b.dataset.topo === state.topology));
}

function selectRecord(i) {
  state.active = i;
  state.selectedOrfId = null;
  state.alignResult = null;
  const rec = activeRec();
  state.topology = rec.topology === 'circular' ? 'circular' : 'linear';
  applyTopologyButtons();
  const A = rec.seq.slice(0, 2000);
  $('#alnA').value = A;
  $('#alnB').value = mutate(A, rec.seq.length);
  state.alignAuto = true;
  renderRecordList();
  refreshAll();
}

function setRecords(recs, idx = 0, syncInput = true) {
  state.records = recs;
  state.active = clamp(idx, 0, recs.length - 1);
  state.selectedOrfId = null;
  state.alignResult = null;
  const rec = activeRec();
  if (rec) {
    state.topology = rec.topology === 'circular' ? 'circular' : 'linear';
    applyTopologyButtons();
    if (syncInput) {
      $('#seqInput').value = rec.seq.length > 20000
        ? `>${rec.name} ${rec.desc}\n` + rec.seq.slice(0, 20000) + '\n...'
        : `>${rec.name} ${rec.desc}\n` + rec.seq;
    }
    const A = rec.seq.slice(0, 2000);
    $('#alnA').value = A;
    $('#alnB').value = mutate(A, rec.seq.length);
    state.alignAuto = true;
  }
  renderRecordList();
  refreshAll();
}

function renderRecordList() {
  const ul = $('#recordList');
  $('#recCount').textContent = state.records.length + ' 条';
  ul.innerHTML = state.records.map((r, i) => `
    <li class="rec ${i === state.active ? 'active' : ''}" data-i="${i}">
      <span class="rec-idx">${i + 1}</span>
      <span class="rec-main">
        <div class="rec-name">${r.name}</div>
        <div class="rec-sub">${r.desc || (r.topology === 'circular' ? '环状' : '线性')}</div>
      </span>
      <span class="rec-bp">${fmt(r.seq.length)}</span>
    </li>`).join('');
  $$('#recordList .rec').forEach((li) => {
    li.addEventListener('click', () => selectRecord(Number(li.dataset.i)));
  });
}

function loadFromText(text) {
  const recs = parseFasta(text);
  if (!recs.length) { showToast('没有解析到有效序列'); return; }
  setRecords(recs, 0, false);
  if (recs.length > 1) showToast(`已解析 ${recs.length} 条序列`);
}

function bindEvents() {
  $$('.tab').forEach((t) => {
    t.addEventListener('click', () => {
      $$('.tab').forEach((x) => x.classList.toggle('active', x === t));
      state.tab = t.dataset.tab;
      $$('.pane').forEach((p) => p.classList.toggle('active', p.dataset.pane === state.tab));
      requestAnimationFrame(renderActiveTab);
    });
  });

  $('#btnParse').addEventListener('click', () => loadFromText($('#seqInput').value));

  const fileInput = $('#fileInput');
  const dz = $('#dropZone');
  dz.addEventListener('click', () => fileInput.click());
  const readFile = (f) => {
    const fr = new FileReader();
    fr.onload = () => {
      $('#seqInput').value = String(fr.result);
      loadFromText(String(fr.result));
    };
    fr.readAsText(f);
  };
  fileInput.addEventListener('change', () => { if (fileInput.files[0]) readFile(fileInput.files[0]); });
  ['dragenter', 'dragover'].forEach((ev) => dz.addEventListener(ev, (e) => { e.preventDefault(); dz.classList.add('hot'); }));
  ['dragleave', 'drop'].forEach((ev) => dz.addEventListener(ev, (e) => { e.preventDefault(); dz.classList.remove('hot'); }));
  dz.addEventListener('drop', (e) => {
    const f = e.dataTransfer.files[0];
    if (f) readFile(f);
  });

  const sel = $('#sampleSel');
  SAMPLE_SEQUENCES.forEach((s) => {
    const o = document.createElement('option');
    o.value = s.id;
    o.textContent = `${s.label} · ${s.seq.length} bp`;
    sel.appendChild(o);
  });
  const loadSample = (id) => {
    const s = SAMPLE_SEQUENCES.find((x) => x.id === id);
    if (!s) return;
    sel.value = '';
    setRecords([{ name: s.name, desc: s.desc, seq: s.seq, topology: s.topology }], 0, true);
  };
  sel.addEventListener('change', () => { if (sel.value) loadSample(sel.value); });
  $('#btnDemo').addEventListener('click', () => loadSample('phix174'));

  const minAA = $('#minAA');
  const syncRange = (input, out, suffix) => {
    const v = Number(input.value);
    out.textContent = v + suffix;
    const p = ((v - Number(input.min)) / (Number(input.max) - Number(input.min))) * 100;
    input.style.setProperty('--fill', p + '%');
  };
  syncRange(minAA, $('#minAAOut'), ' aa');
  syncRange($('#gcWin'), $('#gcWinOut'), ' bp');
  let t1 = null;
  minAA.addEventListener('input', () => {
    syncRange(minAA, $('#minAAOut'), ' aa');
    clearTimeout(t1);
    t1 = setTimeout(refreshAll, 140);
  });
  $('#gcWin').addEventListener('input', () => {
    syncRange($('#gcWin'), $('#gcWinOut'), ' bp');
    renderActiveTab();
  });
  $('#primeOnly').addEventListener('change', () => {
    const rec = activeRec();
    if (rec && rec._siteCache) rec._siteCache.clear();
    refreshAll();
  });

  $('#enzymeSet').addEventListener('change', refreshAll);
  $('#labelMode').addEventListener('change', () => { if (state.tab === 'enzyme') renderEnzyme(); });
  $$('#topoSeg button').forEach((b) => {
    b.addEventListener('click', () => {
      state.topology = b.dataset.topo;
      applyTopologyButtons();
      const rec = activeRec();
      if (rec && rec._siteCache) rec._siteCache.clear();
      refreshAll();
    });
  });

  $$('#orfStrand button').forEach((b) => {
    b.addEventListener('click', () => {
      state.strandFilter = b.dataset.strand;
      $$('#orfStrand button').forEach((x) => x.classList.toggle('active', x === b));
      renderOrfMap();
    });
  });

  $('#orfPicker').addEventListener('change', renderCodon);

  $$('#alnMode button').forEach((b) => {
    b.addEventListener('click', () => {
      $$('#alnMode button').forEach((x) => x.classList.toggle('active', x === b));
    });
  });
  ['#alnA', '#alnB'].forEach((s) => $(s).addEventListener('input', () => { state.alignAuto = false; }));
  $('#btnAlign').addEventListener('click', renderAlign);
  $('#btnMutate').addEventListener('click', () => {
    const A = $('#alnA').value.replace(/[^A-Za-z]/g, '');
    if (!A) { showToast('请先填写序列 A'); return; }
    $('#alnB').value = mutate(A, Date.now() & 0xffff);
    renderAlign();
  });

  const mapCanvas = $('#mapCanvas');
  mapCanvas.addEventListener('mousemove', (e) => {
    const r = mapCanvas.getBoundingClientRect();
    const x = e.clientX - r.left, y = e.clientY - r.top;
    let best = null, bd = Infinity;
    for (const h of state.mapHits) {
      const d = (h.x - x) ** 2 + (h.y - y) ** 2;
      if (d < h.r * h.r && d < bd) { bd = d; best = h; }
    }
    if (best) {
      const s = best.site;
      const end = enzymeEnd(s.enz);
      showTip(`<b>${s.enz.n}</b> <span class="tt-mono">${s.enz.s}</span><br>
        识别位点 <span class="tt-mono">${fmt(s.pos)}</span><br>
        切点 <span class="tt-mono">${fmt(s.cut)}</span> · ${end.t}`, e.clientX, e.clientY);
    } else hideTip();
  });
  mapCanvas.addEventListener('mouseleave', hideTip);

  let rt = null;
  window.addEventListener('resize', () => {
    clearTimeout(rt);
    rt = setTimeout(renderActiveTab, 160);
  });

  document.addEventListener('keydown', (e) => {
    if (e.target.matches('input, textarea, select')) return;
    const map = { 1: 'overview', 2: 'orf', 3: 'enzyme', 4: 'codon', 5: 'align' };
    if (map[e.key]) {
      const tab = $$('.tab').find((t) => t.dataset.tab === map[e.key]);
      if (tab) tab.click();
    }
  });
}

/* ---------------- 启动 ---------------- */
(function boot() {
  bindEvents();
  const first = SAMPLE_SEQUENCES[0];
  if (first) {
    setRecords([{ name: first.name, desc: first.desc, seq: first.seq, topology: first.topology }], 0, true);
  }
  // 深链：#enzyme 直接打开对应标签
  const hash = decodeURIComponent(location.hash.replace('#', ''));
  const tabs = ['overview', 'orf', 'enzyme', 'codon', 'align'];
  if (tabs.includes(hash)) {
    const tab = $$('.tab').find((t) => t.dataset.tab === hash);
    if (tab) tab.click();
  }
  window.addEventListener('hashchange', () => {
    const h = location.hash.replace('#', '');
    const tab = $$('.tab').find((t) => t.dataset.tab === h);
    if (tab) tab.click();
  });
})();
