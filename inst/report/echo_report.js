/* ECHO interactive report -- client-side app.
   Reads the JSON blob in #echo-data, builds one lazily-rendered panel per tab and
   draws every figure with the inlined plotly.js. No server and no network access
   are needed. Architecture follows the SPARKS unified report. */
(function () {
  'use strict';

  // ------------------------------------------------------------------ data & utils
  var DATA = JSON.parse(document.getElementById('echo-data').textContent);
  var META = DATA.meta, TH = DATA.thresholds || {};
  var $ = function (s, r) { return (r || document).querySelector(s); };
  var $$ = function (s, r) { return Array.prototype.slice.call((r || document).querySelectorAll(s)); };
  var HAS_PLOTLY = typeof window.Plotly !== 'undefined';

  function h(tag, attrs) {
    var el = document.createElement(tag);
    if (attrs) Object.keys(attrs).forEach(function (k) {
      var v = attrs[k];
      if (v === null || v === undefined || v === false) return;
      if (k === 'class') el.className = v;
      else if (k === 'text') el.textContent = v;
      else if (k.slice(0, 2) === 'on' && typeof v === 'function') el.addEventListener(k.slice(2), v);
      else if (k === 'style' && typeof v === 'object') Object.assign(el.style, v);
      else el.setAttribute(k, v === true ? '' : v);
    });
    for (var i = 2; i < arguments.length; i++) append(el, arguments[i]);
    return el;
  }
  function append(el, kid) {
    if (kid === null || kid === undefined || kid === false) return;
    if (Array.isArray(kid)) kid.forEach(function (k) { append(el, k); });
    else if (kid.nodeType) el.appendChild(kid);
    else el.appendChild(document.createTextNode(String(kid)));
  }
  function clear(el) { while (el.firstChild) el.removeChild(el.firstChild); return el; }
  function range(n) { var a = new Array(n); for (var i = 0; i < n; i++) a[i] = i; return a; }
  function uniq(a) { return a.filter(function (v, i) { return a.indexOf(v) === i; }); }
  function sum(a) { return a.reduce(function (s, v) { return s + (isFinite(v) ? v : 0); }, 0); }
  function fmtInt(n) { return Number(n).toLocaleString('en-US'); }
  function fmt(v) {
    if (v === null || v === undefined || (typeof v === 'number' && !isFinite(v))) return '';
    if (typeof v !== 'number') return String(v);
    if (Number.isInteger(v)) return String(v);
    if (v !== 0 && Math.abs(v) < 1e-3) return v.toExponential(2);
    return String(Number(v.toPrecision(4)));
  }
  function download(name, text, type) {
    var blob = new Blob([text], { type: type || 'text/csv;charset=utf-8' });
    var a = h('a', { href: URL.createObjectURL(blob), download: name });
    document.body.appendChild(a); a.click();
    setTimeout(function () { URL.revokeObjectURL(a.href); a.remove(); }, 0);
  }
  function csvField(v) {
    var s = v === null || v === undefined ? '' : String(v);
    return /[",\n\r\t]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
  }
  function toCSV(cols, rowFn, n) {
    var lines = [cols.map(csvField).join(',')];
    for (var i = 0; i < n; i++) lines.push(cols.map(function (c) { return csvField(rowFn(c, i)); }).join(','));
    return lines.join('\n') + '\n';
  }
  function store(key, val) {
    try { if (val === undefined) return window.localStorage.getItem(key); window.localStorage.setItem(key, val); } catch (e) { /* storage may be unavailable */ }
    return null;
  }
  function T(tbl) { return tbl && tbl.n ? tbl : { cols: (tbl && tbl.cols) || [], n: 0, data: (tbl && tbl.data) || {} }; }
  function col(tbl, c) { return (tbl && tbl.data && tbl.data[c]) || []; }
  function plural(n, w, ws) { return n + ' ' + (n === 1 ? w : (ws || w + 's')); }

  // ------------------------------------------------------------------ theme & palette
  var CHROME = {
    light: { surface: '#fcfcfb', text: '#0b0b0b', secondary: '#52514e', muted: '#898781', grid: '#e1e0d9', axis: '#c3c2b7', band: '#f0efec' },
    dark: { surface: '#1a1a19', text: '#f0efec', secondary: '#c3c2b7', muted: '#898781', grid: '#2c2c2a', axis: '#383835', band: '#232321' }
  };
  var CAT = {
    light: ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948'],
    dark: ['#3987e5', '#d95926', '#199e70', '#c98500', '#d55181', '#008300', '#9085e9', '#e66767']
  };
  // semantic colours (never the only cue: shapes / labels carry the same meaning)
  var SEM = {
    light: { del: '#e34948', dup: '#2a78d6', high: '#006300', medium: '#eda100', low: '#898781', bad: '#c42b2b', warn: '#8a5a00', ok: '#1baf7a', ref: '#a9a8a0', test: '#2a78d6', hit: '#e34948' },
    dark: { del: '#e66767', dup: '#3987e5', high: '#0ca30c', medium: '#fab219', low: '#898781', bad: '#f08080', warn: '#fab219', ok: '#199e70', ref: '#6f6e68', test: '#3987e5', hit: '#e66767' }
  };
  var SEQ = ['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b'];
  var FONT = "system-ui, -apple-system, 'Segoe UI', sans-serif";

  var themePref = store('echo-theme') || 'auto';
  function resolved() {
    if (themePref === 'light' || themePref === 'dark') return themePref;
    return window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
  }
  function applyThemeAttr() {
    if (themePref === 'auto') document.documentElement.removeAttribute('data-theme');
    else document.documentElement.setAttribute('data-theme', themePref);
  }
  function tint(hex, amt, toward) {
    if (!amt) return hex;
    var n = function (s) { return [1, 3, 5].map(function (i) { return parseInt(s.slice(i, i + 2), 16); }); };
    var a = n(hex), b = n(toward);
    return '#' + a.map(function (v, i) { return ('0' + Math.round(v + (b[i] - v) * amt).toString(16)).slice(-2); }).join('');
  }
  function catColors(n, theme) {
    var base = CAT[theme], out = [], tints = [0, 0.38, 0.62];
    for (var i = 0; i < n; i++) out.push(tint(base[i % 8], tints[Math.min(Math.floor(i / 8), 2)], theme === 'light' ? '#ffffff' : '#000000'));
    return out;
  }

  // ------------------------------------------------------------------ plotly helpers
  var registry = [];
  var CONFIG = { displaylogo: false, responsive: true, modeBarButtonsToRemove: ['autoScale2d'] };

  function baseLayout(t, o) {
    o = o || {};
    var c = CHROME[t];
    var ax = function (label, extra) {
      return Object.assign({
        title: { text: label || '', font: { color: c.secondary } },
        gridcolor: c.grid, zerolinecolor: c.axis, linecolor: c.axis,
        tickfont: { color: c.muted }, automargin: true
      }, extra || {});
    };
    return {
      paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
      font: { family: FONT, size: 12, color: c.secondary },
      xaxis: ax(o.xlab, o.xaxis), yaxis: ax(o.ylab, o.yaxis),
      showlegend: !!o.legend,
      legend: Object.assign({ orientation: 'h', y: 1.12, x: 0, font: { color: c.secondary } }, o.legendOpts || {}),
      margin: o.margin || { l: 56, r: 16, t: 16, b: 48 },
      hoverlabel: { font: { family: FONT } },
      annotations: o.annotations, shapes: o.shapes, barmode: o.barmode, dragmode: o.dragmode, bargap: o.bargap,
      hovermode: o.hovermode
    };
  }
  // build(theme) -> {data, layout, config}
  function plot(div, build) {
    if (!HAS_PLOTLY) { showNote(div, 'plotly.js could not be loaded (offline use needs the report to be generated with report_plotly: inline).'); return div; }
    if (div._note) { clear(div); div._note = false; }
    div._render = function () {
      var s = build(resolved());
      var cfg = Object.assign({}, CONFIG, s.config || {});
      s.layout.autosize = true;
      window.Plotly.react(div, s.data, s.layout, cfg);
      div._dirty = false;
    };
    if (registry.indexOf(div) < 0) registry.push(div);
    div._render();
    return div;
  }
  function rerenderAll() {
    registry = registry.filter(function (d) { return document.body.contains(d) && d._render; });
    registry.forEach(function (d) { if (d.offsetParent !== null) d._render(); else d._dirty = true; });
  }
  function showNote(div, text) {
    if (HAS_PLOTLY && div.data) window.Plotly.purge(div);
    clear(div); div._render = null; div._note = true; div._dirty = false;
    div.appendChild(emptyNote(text));
  }
  function refreshVisible(root) {
    $$('.plot', root).forEach(function (d) {
      if (!d._render) return;
      if (d._dirty) d._render(); else if (HAS_PLOTLY && d.data) window.Plotly.Plots.resize(d);
    });
  }
  function plotBox(cls) { return h('div', { class: 'plot ' + (cls || '') }); }
  function emptyNote(text) { return h('div', { class: 'empty', text: text }); }

  // ------------------------------------------------------------------ tables
  function colType(tb, c) {
    var arr = tb.data[c] || [];
    for (var i = 0; i < arr.length; i++) { if (arr[i] !== null && arr[i] !== undefined) return typeof arr[i]; }
    return 'string';
  }
  // tb = {cols, data: {col: [..]}, n}. opts: filename, pageSize, render: {col: fn(value,i)->node}, sort: [col,'asc'|'desc'],
  // onRow(i), wrap: [cols], note
  function tableView(host, tb, opts) {
    opts = opts || {};
    var n = tb.n, pageSize = opts.pageSize || 25, page = 0, filter = '', pred = null, selIdx = -1;
    var sortCol = opts.sort ? opts.sort[0] : null, sortDir = opts.sort ? opts.sort[1] : 'asc';
    var types = {}; tb.cols.forEach(function (c) { types[c] = colType(tb, c); });
    var idx = range(n);
    var search = h('input', { type: 'search', placeholder: 'Filter rows…', 'aria-label': 'Filter rows' });
    var count = h('span', { class: 'muted small' });
    var pager = h('span', { class: 'pager' });
    var tbody = h('tbody'), thead = h('thead');
    var table = h('table', { class: 'tbl' }, thead, tbody);
    var wrap = h('div', { class: 'tbl-wrap' }, table);
    var dl = h('button', { class: 'btn', type: 'button', text: 'Download CSV', onclick: function () {
      download((opts.filename || 'table') + '.csv', toCSV(tb.cols, function (c, k) { return tb.data[c][idx[k]]; }, idx.length));
    } });
    host.appendChild(h('div', { class: 'tbl-bar' }, search, count, dl, pager));
    host.appendChild(wrap);
    if (opts.note) host.appendChild(h('div', { class: 'note', text: opts.note }));

    function recompute() {
      var f = filter.toLowerCase();
      idx = [];
      for (var i = 0; i < n; i++) {
        if (pred && !pred(i)) continue;
        if (f) {
          var hit = false;
          for (var j = 0; j < tb.cols.length && !hit; j++) {
            var v = tb.data[tb.cols[j]][i];
            if (v !== null && v !== undefined && String(v).toLowerCase().indexOf(f) >= 0) hit = true;
          }
          if (!hit) continue;
        }
        idx.push(i);
      }
      if (sortCol) {
        var arr = tb.data[sortCol], dir = sortDir === 'asc' ? 1 : -1, num = types[sortCol] === 'number';
        idx.sort(function (a, b) {
          var x = arr[a], y = arr[b];
          if (x === null || x === undefined) return 1;
          if (y === null || y === undefined) return -1;
          return (num ? x - y : String(x).localeCompare(String(y), undefined, { numeric: true })) * dir;
        });
      }
      page = Math.min(page, Math.max(0, Math.ceil(idx.length / pageSize) - 1));
      draw();
    }
    function draw() {
      var tr = h('tr');
      tb.cols.forEach(function (c) {
        tr.appendChild(h('th', {
          scope: 'col', tabindex: '0', text: c, 'aria-sort': sortCol === c ? (sortDir === 'asc' ? 'ascending' : 'descending') : 'none',
          onclick: function () { sortDir = sortCol === c && sortDir === 'asc' ? 'desc' : 'asc'; sortCol = c; recompute(); },
          onkeydown: function (e) { if (e.key === 'Enter') this.click(); }
        }));
      });
      clear(thead).appendChild(tr);
      clear(tbody);
      var start = page * pageSize, end = Math.min(idx.length, start + pageSize);
      var wrapCols = opts.wrap || [];
      for (var k = start; k < end; k++) {
        var i = idx[k], row = h('tr', opts.onRow ? { class: 'pick' + (i === selIdx ? ' sel' : ''), tabindex: '0' } : null);
        if (opts.onRow) (function (ii, rr) {
          rr.addEventListener('click', function () { selIdx = ii; $$('tr.sel', tbody).forEach(function (e) { e.classList.remove('sel'); }); rr.classList.add('sel'); opts.onRow(ii); });
          rr.addEventListener('keydown', function (e) { if (e.key === 'Enter') rr.click(); });
        })(i, row);
        tb.cols.forEach(function (c) {
          var v = tb.data[c][i], r = opts.render && opts.render[c];
          row.appendChild(h('td', { class: (types[c] === 'number' ? 'num' : '') + (wrapCols.indexOf(c) >= 0 ? ' wrap' : '') }, r ? r(v, i) : fmt(v)));
        });
        tbody.appendChild(row);
      }
      count.textContent = fmtInt(idx.length) + (idx.length === n ? '' : ' of ' + fmtInt(n)) + ' rows';
      var pages = Math.max(1, Math.ceil(idx.length / pageSize));
      clear(pager);
      pager.appendChild(h('button', { class: 'btn', type: 'button', text: '‹', 'aria-label': 'Previous page', disabled: page === 0, onclick: function () { page--; draw(); } }));
      pager.appendChild(h('span', { class: 'small', text: (page + 1) + ' / ' + pages }));
      pager.appendChild(h('button', { class: 'btn', type: 'button', text: '›', 'aria-label': 'Next page', disabled: page >= pages - 1, onclick: function () { page++; draw(); } }));
    }
    var timer = null;
    search.addEventListener('input', function () {
      clearTimeout(timer);
      timer = setTimeout(function () { filter = search.value; page = 0; recompute(); }, 120);
    });
    recompute();
    return {
      setPredicate: function (fn) { pred = fn; page = 0; recompute(); },
      refresh: recompute,
      visible: function () { return idx.slice(); },
      select: function (i) { selIdx = i; draw(); }
    };
  }

  // ------------------------------------------------------------------ global state (sample selector)
  var state = { sample: META.focus || '', listeners: [] };
  function sel() { return state.sample; }
  function setSample(s, quiet) {
    state.sample = s || '';
    var el = $('#sample-sel'); if (el) el.value = state.sample;
    if (!quiet) state.listeners.forEach(function (cb) { cb(); });
  }
  function onSample(cb) { state.listeners.push(cb); }
  function sampleRow(S, name) { var i = col(S, 'Sample').indexOf(name); return i < 0 ? null : i; }
  function chip(text, cls) { return h('span', { class: 'chip ' + (cls || '') , text: text }); }
  function confChip(v) { return v ? chip(v, v) : ''; }
  function sevChip(v) { return v ? chip(v, v) : ''; }

  // ------------------------------------------------------------------ alerts
  var SEV_ICON = { error: '✖', warning: '▲', info: 'ℹ' };
  var SEV_ORDER = { error: 0, warning: 1, info: 2 };
  function activeAlerts(kinds, minSev) {
    var s = sel();
    return (DATA.alerts || []).filter(function (a) {
      if (kinds && kinds.indexOf(a.kind) < 0) return false;
      if (minSev && SEV_ORDER[a.severity] > SEV_ORDER[minSev]) return false;
      if (s && a.samples && a.samples.length && a.samples.indexOf(s) < 0) return false;
      return true;
    }).sort(function (a, b) { return SEV_ORDER[a.severity] - SEV_ORDER[b.severity]; });
  }
  function alertsBox(kinds, minSev, withLinks) {
    var list = activeAlerts(kinds, minSev);
    if (!list.length) return null;
    return h('div', { class: 'alerts', role: 'region', 'aria-label': 'Warnings' }, list.map(function (a) {
      return h('div', { class: 'alert ' + a.severity, role: a.severity === 'info' ? 'note' : 'alert' },
        h('span', { class: 'ic', 'aria-hidden': 'true', text: SEV_ICON[a.severity] }),
        h('div', { class: 'body' },
          h('div', { class: 'ttl', text: a.title }),
          h('div', { text: a.text }),
          (withLinks !== false && a.tab && a.tab !== active) ? h('a', { href: '#tab=' + a.tab, text: 'Details →', onclick: function (e) { e.preventDefault(); nav(a.tab); } }) : null));
    }));
  }

  // ------------------------------------------------------------------ tabs
  var TABS = [];
  function tab(id, label, available, build) { TABS.push({ id: id, label: label, available: available, build: build, built: false }); }
  function card(title, hint, body) {
    return h('section', { class: 'card' }, h('header', null, h('h3', { text: title }), hint ? h('span', { class: 'muted small hint', text: hint }) : null), body);
  }
  function tile(label, value, cls, sub) {
    return h('div', { class: 'tile ' + (cls || '') }, h('div', { class: 'v', text: String(value) }), h('div', { class: 'l', text: label }), sub ? h('div', { class: 'note', text: sub }) : null);
  }
  function nav(id) { activate(id, true); }
  function filterBySample(tbl, key) {
    return function (i) { var s = sel(); if (!s) return true; var v = tbl.data[key || 'Sample'][i]; return v === s || v === 'All'; };
  }
  function sampleOnly(tbl) { return function (i) { var s = sel(); return !s || tbl.data.Sample[i] === s; }; }
  function n_calls() { var C = DATA.cnvs, s = sel(); return s ? col(C, 'Sample').filter(function (x) { return x === s; }).length : C.n; }
  function shortList(a, k) { return a.length <= k ? a.join(', ') : a.slice(0, k).join(', ') + ' … (+' + (a.length - k) + ')'; }

  // ---- Overview
  tab('overview', 'Overview', function () { return true; }, function (root) {
    var S = DATA.samples, C = DATA.cnvs, s = sel();
    var ab = alertsBox(null, 'warning'); if (ab) root.appendChild(ab);

    var analysed = col(S, 'Analysed').filter(function (v) { return v === 'Yes'; }).length;
    var rows = range(S.n).filter(function (i) { return !s || S.data.Sample[i] === s; });
    var flagged = rows.filter(function (i) { return S.data.Status[i] === 'Warning' || S.data.Status[i] === 'Fail'; }).length;
    var calls = range(C.n).filter(function (i) { return !s || C.data.Sample[i] === s; });
    var high = calls.filter(function (i) { return C.data.Confidence[i] === 'HIGH'; }).length;
    var scopeTxt = META.modechrom_label;
    var tiles = [
      tile(s ? 'Sample' : 'Samples analysed', s ? s : analysed + ' / ' + S.n, ''),
      tile('Targets (exons)', fmtInt(META.n_targets), ''),
      tile('CNV calls', calls.length, calls.length ? '' : 'ok'),
      tile('HIGH confidence', high, ''),
      tile(s ? 'QC status' : 'Samples with QC flags', s ? (rows.length ? S.data.Status[rows[0]] : '') : flagged, flagged ? 'warn' : 'ok'),
      tile('Chromosomes analysed', scopeTxt, META.scope_incomplete ? 'warn' : 'ok', META.scope_incomplete ? 'The remaining chromosome classes are NOT analysed in this run' : null)
    ];
    root.appendChild(h('div', { class: 'tiles', style: { marginBottom: '12px' } }, tiles));

    var left = plotBox('short'), right = plotBox('short');
    var g = h('div', { class: 'grid g2' }, card('CNV calls per sample', 'click a bar to filter the report to that sample', left),
      card('CNV landscape', 'log2(observed / expected) along the targets; marker size = number of exons', right));
    root.appendChild(g);
    setTimeout(function () { callsPerSample(left); landscape(right); }, 0);

    var attention = rows.filter(function (i) { return S.data.Status[i] !== 'OK'; });
    if (!s && attention.length) {
      var host = h('div', { class: 'pill-list' });
      attention.slice(0, 40).forEach(function (i) {
        host.appendChild(h('button', { class: 'btn', type: 'button', title: S.data['QC flags'][i] || '', onclick: function () { setSample(S.data.Sample[i]); nav('samples'); } },
          S.data.Sample[i] + ' ', chip(S.data.Status[i], S.data.Status[i] === 'Fail' ? 'fail' : S.data.Status[i] === 'Warning' ? 'warn' : 'info')));
      });
      root.appendChild(h('div', { style: { marginTop: '12px' } }, card('Samples needing attention', attention.length + ' flagged or excluded', host)));
    }
  });

  function callsPerSample(div) {
    var S = DATA.samples, C = DATA.cnvs;
    if (!S.n) { showNote(div, 'No samples.'); return; }
    var names = S.data.Sample.slice(), s = sel();
    if (names.length > 60) names = names.filter(function (nm) { return col(C, 'Sample').indexOf(nm) >= 0; });
    var levels = ['HIGH', 'MEDIUM', 'LOW'];
    plot(div, function (t) {
      var sem = SEM[t], colr = { HIGH: sem.high, MEDIUM: sem.medium, LOW: sem.low };
      var data = levels.map(function (lv) {
        return { type: 'bar', name: lv, x: names, y: names.map(function (nm) {
          var k = 0; for (var i = 0; i < C.n; i++) if (C.data.Sample[i] === nm && C.data.Confidence[i] === lv) k++; return k; }),
          marker: { color: colr[lv], opacity: names.map(function (nm) { return !s || nm === s ? 1 : 0.3; }) },
          hovertemplate: '%{x}<br>' + lv + ': %{y}<extra></extra>' };
      });
      var L = baseLayout(t, { ylab: 'CNV calls', legend: true, barmode: 'stack', xaxis: { tickangle: -45, type: 'category' }, margin: { l: 56, r: 16, t: 28, b: 90 } });
      L.yaxis.rangemode = 'tozero';
      var top = Math.max.apply(null, names.map(function (nm) { var k = 0; for (var i = 0; i < C.n; i++) if (C.data.Sample[i] === nm) k++; return k; }));
      if (top <= 6) { L.yaxis.dtick = 1; L.yaxis.range = [0, Math.max(top, 1) + 0.3]; }
      return { data: data, layout: L };
    });
    if (div.on) div.on('plotly_click', function (ev) { var p = ev.points && ev.points[0]; if (p) { setSample(sel() === p.x ? '' : p.x); } });
  }

  function landscape(div) {
    var C = DATA.cnvs, s = sel(), CH = DATA.chr;
    var idx = range(C.n).filter(function (i) { return !s || C.data.Sample[i] === s; });
    if (!idx.length) { showNote(div, 'No CNV calls to show.'); return; }
    plot(div, function (t) {
      var sem = SEM[t], c = CHROME[t];
      var shapes = [], tickv = [], tickt = [];
      range(CH.n).forEach(function (k) {
        var a = CH.data.first[k] - 0.5, b = CH.data.last[k] + 0.5;
        if (k % 2 === 0) shapes.push({ type: 'rect', xref: 'x', yref: 'paper', x0: a, x1: b, y0: 0, y1: 1, fillcolor: c.band, line: { width: 0 }, layer: 'below' });
        tickv.push((a + b) / 2); tickt.push(CH.data.chromosome[k].replace(/^chr/, ''));
      });
      var data = ['DEL', 'DUP'].map(function (ty) {
        var ii = idx.filter(function (i) { return C.data.Type[i] === ty; });
        return { type: 'scatter', mode: 'markers', name: ty === 'DEL' ? 'Deletion' : 'Duplication',
          x: ii.map(function (i) { return C.data.gidx[i]; }),
          y: ii.map(function (i) { return C.data.Ratio[i] > 0 ? Math.log2(C.data.Ratio[i]) : null; }),
          text: ii.map(function (i) { return C.data.Sample[i] + '<br>' + C.data.Gene[i] + ' (' + C.data.Genomic[i] + ')<br>' + C.data.Confidence[i] + ' · ratio ' + fmt(C.data.Ratio[i]); }),
          hoverinfo: 'text',
          marker: { color: ty === 'DEL' ? sem.del : sem.dup, symbol: ty === 'DEL' ? 'triangle-down' : 'triangle-up',
            size: ii.map(function (i) { return Math.min(22, 8 + 2 * Math.sqrt(C.data.Nexons[i] || 1)); }), line: { width: 1, color: c.surface } } };
      });
      var L = baseLayout(t, { ylab: 'log2(obs / exp)', legend: true, shapes: shapes, xaxis: { tickmode: 'array', tickvals: tickv, ticktext: tickt, tickangle: 0, showgrid: false } });
      L.yaxis.zeroline = true;
      return { data: data, layout: L };
    });
  }

  // ---- QC & warnings
  tab('qc', 'QC & warnings', function () { return true; }, function (root) {
    var S = DATA.samples, Q = DATA.qc, CH = DATA.chr;
    var ab = alertsBox(null, 'info', false); if (ab) root.appendChild(ab);
    // chromosome scope
    var sc = DATA.scope, scopeTiles = range(sc.n).filter(function (i) { return sc.data.class[i] !== 'other' || sc.data.present[i]; }).map(function (i) {
      var present = sc.data.present[i], analysed = sc.data.analysed[i];
      var status = !present ? 'not in BED' : analysed ? 'analysed' : 'NOT analysed';
      var cls = !present ? 'warn' : analysed ? 'ok' : 'bad';
      if (sc.data.class[i] === 'other') cls = '';
      return tile(sc.data.label[i], fmtInt(sc.data.targets[i]) + ' targets', cls, status);
    });
    root.appendChild(card('Chromosome scope of this run', 'modechrom = ' + META.modechrom + ' (' + META.modechrom_label + ')',
      h('div', null,
        h('p', { class: 'note', text: 'ECHO calls CNVs on one chromosome class per run: autosomes (A), chrX (X) or chrY (Y). Classes that are present in the BED but not analysed here need their own run (change modechrom).' }),
        h('div', { class: 'scope' }, scopeTiles))));

    var m1 = plotBox('short'), m2 = plotBox('short'), m3 = plotBox('short'), m4 = plotBox('short');
    root.appendChild(h('div', { class: 'grid g2', style: { marginTop: '12px' } },
      card('Targets per chromosome', 'colour = analysed in this run or not', m1),
      card('Median depth per sample', 'threshold min_cov = ' + fmt(TH.min_cov), m2)));
    root.appendChild(h('div', { class: 'grid g2', style: { marginTop: '12px' } },
      card('Total reads per sample', 'threshold min_total_reads = ' + fmtInt(TH.min_total_reads), m3),
      card('Best correlation with another sample', 'threshold min_corr = ' + fmt(TH.min_corr), m4)));
    setTimeout(function () {
      chromBars(m1);
      sampleMetric(m2, 'Median depth', TH.min_cov, 'Median depth', false);
      sampleMetric(m3, 'Total reads', TH.min_total_reads, 'Total reads', false);
      sampleMetric(m4, 'Max correlation', TH.min_corr, 'Spearman r', true);
    }, 0);

    // flags table
    var host = h('div'), sev = '', typ = '';
    var tv = null;
    var bar = h('div', { class: 'toolbar' });
    var segBtns = ['', 'error', 'warning', 'info'].map(function (v) {
      return h('button', { type: 'button', 'aria-pressed': v === sev ? 'true' : 'false', text: v ? v[0].toUpperCase() + v.slice(1) : 'All', onclick: function () {
        sev = v; segBtns.forEach(function (b, k) { b.setAttribute('aria-pressed', ['', 'error', 'warning', 'info'][k] === v ? 'true' : 'false'); }); apply(); } });
    });
    var types = uniq(col(Q, 'Type')).sort();
    var typeSel = h('select', { 'aria-label': 'Flag type' }, h('option', { value: '', text: 'All types' }), types.map(function (x) { return h('option', { value: x, text: x }); }));
    typeSel.addEventListener('change', function () { typ = typeSel.value; apply(); });
    bar.appendChild(h('label', { class: 'field' }, 'Severity', h('div', { class: 'seg', role: 'group' }, segBtns)));
    bar.appendChild(h('label', { class: 'field' }, 'Type', typeSel));
    host.appendChild(bar);
    var tblHost = h('div'); host.appendChild(tblHost);
    function apply() {
      if (!tv) return;
      var bySample = filterBySample(Q);
      tv.setPredicate(function (i) { return bySample(i) && (!sev || Q.data.Severity[i] === sev) && (!typ || Q.data.Type[i] === typ); });
    }
    if (!Q.n) {
      root.appendChild(h('div', { style: { marginTop: '12px' } }, card('QC flags', null, emptyNote(DATA.meta.qc_available ? 'No QC flag was raised.' : 'No QC metrics file was found for this run.'))));
    } else {
      tv = tableView(tblHost, Q, { filename: 'echo_qc_flags', pageSize: 20, wrap: ['Details'],
        sort: ['Severity', 'asc'], render: { Severity: function (v) { return sevChip(v); }, Details: function (v) { return v; } } });
      apply();
      root.appendChild(h('div', { style: { marginTop: '12px' } }, card('QC flags', 'Whole sample · Whole exon · Exon variability · Cohort size · Missing Chromosome · Analysis scope', host)));
    }
  });

  function chromBars(div) {
    var CH = DATA.chr;
    if (!CH.n) { showNote(div, 'No targets.'); return; }
    plot(div, function (t) {
      var sem = SEM[t];
      var colors = range(CH.n).map(function (i) { return CH.data.analysed[i] ? sem.ok : sem.low; });
      var data = [{ type: 'bar', x: CH.data.chromosome, y: CH.data.targets, marker: { color: colors },
        customdata: range(CH.n).map(function (i) { return CH.data.analysed[i] ? 'analysed' : 'NOT analysed'; }),
        hovertemplate: '%{x}: %{y} targets (%{customdata})<extra></extra>' }];
      var L = baseLayout(t, { ylab: 'Targets', xaxis: { type: 'category', tickangle: -45 }, margin: { l: 56, r: 16, t: 16, b: 60 } });
      return { data: data, layout: L };
    });
  }
  function sampleMetric(div, colName, threshold, ylab, higherOk) {
    var S = DATA.samples;
    var rows = range(S.n).filter(function (i) { return S.data[colName][i] !== null && S.data[colName][i] !== undefined; });
    if (!rows.length) { showNote(div, 'Not available (needs at least 2 samples).'); return; }
    plot(div, function (t) {
      var sem = SEM[t];
      var vals = rows.map(function (i) { return S.data[colName][i]; });
      var bad = vals.map(function (v) { return threshold !== null && threshold !== undefined && v < threshold; });
      var s = sel();
      var data = [{ type: 'bar', x: rows.map(function (i) { return S.data.Sample[i]; }), y: vals,
        marker: { color: bad.map(function (b) { return b ? sem.bad : sem.dup; }), line: { width: rows.map(function (i) { return s && S.data.Sample[i] === s ? 2 : 0; }), color: CHROME[t].text } },
        hovertemplate: '%{x}<br>' + ylab + ': %{y}<extra></extra>' }];
      var shapes = (threshold !== null && threshold !== undefined) ? [{ type: 'line', xref: 'paper', x0: 0, x1: 1, yref: 'y', y0: threshold, y1: threshold, line: { color: sem.warn, width: 2, dash: 'dash' } }] : [];
      var L = baseLayout(t, { ylab: ylab, shapes: shapes, xaxis: { type: 'category', tickangle: -45 }, margin: { l: 64, r: 16, t: 16, b: 90 } });
      if (higherOk) { var lo = Math.min.apply(null, vals.concat([threshold])); L.yaxis.range = [Math.max(0, lo - 0.02), 1.0]; }
      return { data: data, layout: L };
    });
  }

  // ---- Samples
  tab('samples', 'Samples', function () { return DATA.samples.n > 0; }, function (root) {
    var S = DATA.samples;
    var ab = alertsBox(['sample', 'cohort'], 'info', false); if (ab) root.appendChild(ab);
    var cols = ['Sample', 'Status', 'Analysed'].concat(S.cols.filter(function (c) { return ['Sample', 'Status', 'Analysed', 'QC flags', 'QC severity'].indexOf(c) < 0; })).concat(['QC flags']);
    var disp = { cols: cols, n: S.n, data: S.data };
    var detail = h('div', { style: { marginTop: '12px' } });
    var host = h('div');
    var tv = tableView(host, disp, { filename: 'echo_samples', pageSize: 25, wrap: ['QC flags'], sort: ['Status', 'asc'],
      render: {
        Status: function (v) { return chip(v, v === 'Fail' ? 'fail' : v === 'Warning' ? 'warn' : v === 'Excluded' ? 'info' : 'ok'); },
        Sample: function (v) { return h('a', { href: '#', text: v, onclick: function (e) { e.preventDefault(); setSample(v); } }); }
      },
      onRow: function (i) { showSample(detail, i); } });
    tv.setPredicate(sampleOnly(S));
    root.appendChild(card('Samples', 'click a row for details; the name filters the whole report', host));
    root.appendChild(detail);
    if (sel()) { var k = sampleRow(S, sel()); if (k !== null) showSample(detail, k); }
  });
  function showSample(host, i) {
    var S = DATA.samples, nm = S.data.Sample[i];
    clear(host);
    var kv = h('dl', { class: 'kv' });
    S.cols.forEach(function (c) { if (c === 'QC flags') return; var v = S.data[c][i]; if (v === null || v === undefined || v === '') return; kv.appendChild(h('dt', { text: c })); kv.appendChild(h('dd', { text: fmt(v) })); });
    var flags = col(DATA.qc, 'Sample').map(function (x, k) { return k; }).filter(function (k) { return DATA.qc.data.Sample[k] === nm; });
    var list = flags.length ? h('ul', null, flags.map(function (k) { return h('li', null, sevChip(DATA.qc.data.Severity[k]), ' ', DATA.qc.data.Type[k] + ': ' + DATA.qc.data.Details[k]); })) : h('p', { class: 'note', text: 'No sample-level QC flag.' });
    var refs = (DATA.refs && DATA.refs[nm]) || [];
    host.appendChild(card(nm, S.data.Status[i], h('div', { class: 'grid g2' },
      h('div', null, kv),
      h('div', null, h('h3', { text: 'QC flags' }), list,
        h('h3', { text: 'Reference samples (' + refs.length + ')', style: { marginTop: '10px' } }), h('p', { class: 'note', text: refs.length ? refs.join(', ') : 'none' }),
        h('p', null, h('button', { class: 'btn primary', type: 'button', text: 'Show CNV calls', onclick: function () { setSample(nm); nav('cnvs'); } }))))));
  }

  // ---- CNV calls
  tab('cnvs', 'CNV calls', function () { return true; }, function (root) {
    var C = DATA.cnvs;
    var ab = alertsBox(['scope', 'cohort', 'run'], 'warning', true); if (ab) root.appendChild(ab);
    if (!C.n) { root.appendChild(card('CNV calls', null, emptyNote('No CNV call was made in this run.'))); return; }
    var conf = { HIGH: true, MEDIUM: true, LOW: true }, type = '', chr = '', tv = null;
    var tb = h('div', { class: 'toolbar' });
    var confBtns = h('div', { class: 'seg', role: 'group', 'aria-label': 'Confidence' }, ['HIGH', 'MEDIUM', 'LOW'].map(function (k) {
      return h('button', { type: 'button', 'aria-pressed': 'true', text: k, onclick: function () { conf[k] = !conf[k]; this.setAttribute('aria-pressed', conf[k] ? 'true' : 'false'); apply(); } });
    }));
    var typeSel = h('select', { 'aria-label': 'Call type' }, [['', 'All types'], ['DEL', 'Deletions'], ['DUP', 'Duplications']].map(function (o) { return h('option', { value: o[0], text: o[1] }); }));
    typeSel.addEventListener('change', function () { type = typeSel.value; apply(); });
    var chrs = uniq(col(C, 'Chr'));
    var chrSel = h('select', { 'aria-label': 'Chromosome' }, h('option', { value: '', text: 'All chromosomes' }), chrs.map(function (x) { return h('option', { value: x, text: x }); }));
    chrSel.addEventListener('change', function () { chr = chrSel.value; apply(); });
    tb.appendChild(h('label', { class: 'field' }, 'Confidence', confBtns));
    tb.appendChild(h('label', { class: 'field' }, 'Type', typeSel));
    tb.appendChild(h('label', { class: 'field' }, 'Chromosome', chrSel));
    var count = h('span', { class: 'muted small' });
    tb.appendChild(count);
    var host = h('div'), detail = h('div', { style: { marginTop: '12px' } });
    var cols = ['ID', 'Sample', 'Chr', 'Start', 'End', 'Gene', 'Type', 'Exons', 'Nexons', 'Ratio', 'BF', 'Confidence', 'Correlation', 'N.comp', 'Calibration'];
    cols = cols.filter(function (c) { return C.data[c]; });
    var disp = { cols: cols, n: C.n, data: C.data };
    tv = tableView(host, disp, { filename: 'echo_cnv_calls', pageSize: 20, sort: ['ID', 'asc'],
      render: {
        Confidence: function (v) { return confChip(v); },
        Type: function (v) { return chip(v, v); },
        Calibration: function (v) { return v === 'check' ? chip('check', 'warn') : v ? chip(v, 'ok') : ''; },
        Start: function (v) { return fmtInt(v); }, End: function (v) { return fmtInt(v); }
      },
      onRow: function (i) { showCall(detail, i); } });
    function apply() {
      var bySample = sampleOnly(C);
      tv.setPredicate(function (i) { return bySample(i) && conf[C.data.Confidence[i]] !== false && (!type || C.data.Type[i] === type) && (!chr || C.data.Chr[i] === chr); });
      var v = tv.visible().length; count.textContent = plural(v, 'call') + ' shown';
    }
    apply();
    root.appendChild(card('CNV calls', 'click a row to see the coverage evidence', h('div', null, tb, host)));
    root.appendChild(detail);
    var first = tv.visible()[0];
    if (first !== undefined && C.n <= 400) { tv.select(first); showCall(detail, first); }
  });

  function showCall(host, i) {
    var C = DATA.cnvs, id = String(C.data.ID[i]), W = DATA.windows[id];
    clear(host);
    var head = h('div', { class: 'detail-head' },
      h('h2', { text: C.data.Sample[i] + ' · ' + C.data.Gene[i] }),
      chip(C.data.Type[i], C.data.Type[i]), confChip(C.data.Confidence[i]),
      h('span', { class: 'muted', text: C.data.Genomic[i] }));
    var kv = h('dl', { class: 'kv' });
    [['Exons (within gene)', C.data.Exons[i]], ['Number of exons', C.data.Nexons[i]], ['Reads observed', C.data.Observed[i]], ['Reads expected', C.data.Expected[i]],
      ['Observed / expected', C.data.Ratio[i]], ['Bayes factor', C.data.BF[i]], ['Correlation with references', C.data.Correlation[i]],
      ['Reference samples', C.data['N.comp'][i]], ['Comparators', C.data.Comparators[i]]].forEach(function (p) {
      if (p[1] === null || p[1] === undefined || p[1] === '') return; kv.appendChild(h('dt', { text: p[0] })); kv.appendChild(h('dd', { text: fmt(p[1]) }));
    });
    var body = h('div');
    if (!W) {
      body.appendChild(emptyNote(DATA.meta.windows_limited ? 'Coverage profiles were only embedded for the first ' + DATA.meta.windows_limit + ' calls (report_max_profiles).' : 'No coverage window could be computed for this call.'));
    } else {
      if (W.subtitle) body.appendChild(h('div', { class: 'alert warning', role: 'alert' }, h('span', { class: 'ic', text: SEV_ICON.warning }), h('div', { class: 'body' }, h('div', { class: 'ttl', text: 'Background calibration' }), h('div', { text: W.subtitle }))));
      var strip = plotBox('strip'), cov = plotBox('mid'), ci = plotBox('mid'), z = plotBox('mid');
      body.appendChild(strip); body.appendChild(card('Log coverage (scaled references vs test sample)', null, cov));
      body.appendChild(h('div', { style: { height: '12px' } }));
      body.appendChild(card('Observed / expected with 95% beta-binomial interval', null, ci));
      body.appendChild(h('div', { style: { height: '12px' } }));
      body.appendChild(card('Z-score vs references', null, z));
      setTimeout(function () { winPlots(W, strip, cov, ci, z); }, 0);
    }
    host.appendChild(card('Call details', 'CNV.ID ' + id, h('div', null, head, h('div', { class: 'grid g-main' }, h('div', { class: 'stack' }, body), h('div', null, kv)))));
    host.scrollIntoView && host.scrollIntoView({ block: 'nearest', behavior: 'smooth' });
  }

  function winAxis(W, t) {
    var px = W.px, lo = Math.min.apply(null, px) - 0.8, hi = Math.max.apply(null, px) + 0.8;
    return { tickmode: 'array', tickvals: px, ticktext: W.xlab, tickangle: -45, range: [lo, hi], showgrid: false };
  }
  function winPlots(W, strip, cov, ci, z) {
    var px = W.px, aff = W.affected;
    var offRe = null; try { if (META.off_target_pattern) offRe = new RegExp(META.off_target_pattern); } catch (e) { offRe = null; }
    // gene track
    plot(strip, function (t) {
      var cc = catColors(8, t), shapes = [], ann = [], k = 0, start = 0;
      for (var i = 1; i <= px.length; i++) {
        if (i === px.length || W.gene_group[i] !== W.gene_group[start]) {
          var g = W.gene[start] || '', off = offRe && offRe.test(g);
          shapes.push({ type: 'rect', xref: 'x', yref: 'paper', x0: px[start] - 0.5, x1: px[i - 1] + 0.5, y0: 0.05, y1: 0.95, fillcolor: off ? CHROME[t].axis : cc[k % 8], line: { width: 0 }, opacity: 0.9 });
          if (!off && g) ann.push({ x: (px[start] + px[i - 1]) / 2, y: 0.5, xref: 'x', yref: 'paper', text: g, showarrow: false, font: { color: '#ffffff', size: 12 } });
          k++; start = i;
        }
      }
      var L = baseLayout(t, { shapes: shapes, annotations: ann, xaxis: Object.assign(winAxis(W, t), { showticklabels: false }), yaxis: { visible: false }, margin: { l: 56, r: 16, t: 4, b: 4 } });
      L.yaxis = { visible: false, range: [0, 1], fixedrange: true };
      return { data: [{ type: 'scatter', x: [px[0]], y: [0.5], mode: 'markers', marker: { opacity: 0 }, hoverinfo: 'skip' }], layout: L, config: { staticPlot: true } };
    });
    // coverage
    plot(cov, function (t) {
      var sem = SEM[t], data = [];
      var rx = [], ry = [], rt = [];
      W.ref_names.forEach(function (nm, r) { W.ref_cov[r].forEach(function (v, j) { rx.push(px[j]); ry.push(v); rt.push(nm + '<br>exon ' + (W.xlab[j] || '')); }); });
      data.push({ type: 'scatter', mode: 'markers', name: 'Reference samples', x: rx, y: ry, text: rt, hoverinfo: 'text+y', marker: { color: sem.ref, size: 6, opacity: 0.7 } });
      var ni = range(px.length).filter(function (j) { return !aff[j]; }), ai = range(px.length).filter(function (j) { return aff[j]; });
      data.push({ type: 'scatter', mode: 'markers', name: 'Test sample', x: ni.map(function (j) { return px[j]; }), y: ni.map(function (j) { return W.test[j]; }), text: ni.map(function (j) { return 'exon ' + W.xlab[j]; }), hoverinfo: 'text+y', marker: { color: sem.test, size: 8 } });
      data.push({ type: 'scatter', mode: 'markers', name: 'Affected exon(s)', x: ai.map(function (j) { return px[j]; }), y: ai.map(function (j) { return W.test[j]; }), text: ai.map(function (j) { return 'exon ' + W.xlab[j]; }), hoverinfo: 'text+y', marker: { color: sem.hit, size: 11, symbol: 'diamond' } });
      return { data: data, layout: baseLayout(t, { ylab: 'log(coverage)', legend: true, xaxis: winAxis(W, t) }) };
    });
    // observed / expected
    plot(ci, function (t) {
      var sem = SEM[t], data = [], gs = uniq(W.gene_group);
      gs.forEach(function (g, k) {
        var ii = range(px.length).filter(function (j) { return W.gene_group[j] === g; });
        data.push({ type: 'scatter', mode: 'lines', x: ii.map(function (j) { return px[j]; }), y: ii.map(function (j) { return W.lo[j]; }), line: { width: 0 }, hoverinfo: 'skip', showlegend: false });
        data.push({ type: 'scatter', mode: 'lines', x: ii.map(function (j) { return px[j]; }), y: ii.map(function (j) { return W.hi[j]; }), line: { width: 0 }, fill: 'tonexty', fillcolor: CHROME[t].grid, name: '95% interval', showlegend: k === 0, hoverinfo: 'skip' });
      });
      var ni = range(px.length).filter(function (j) { return !aff[j]; }), ai = range(px.length).filter(function (j) { return aff[j]; });
      var tip = function (j) { return 'exon ' + W.xlab[j] + '<br>observed ' + fmt(W.observed[j]) + '<br>expected ' + fmt(W.expected[j]) + '<br>ratio ' + fmt(W.ratio[j]) + ' [' + fmt(W.lo[j]) + ' – ' + fmt(W.hi[j]) + ']'; };
      data.push({ type: 'scatter', mode: 'markers', name: 'Observed', x: ni.map(function (j) { return px[j]; }), y: ni.map(function (j) { return W.ratio[j]; }), text: ni.map(tip), hoverinfo: 'text', marker: { color: sem.test, size: 9 } });
      data.push({ type: 'scatter', mode: 'markers', name: 'Affected', x: ai.map(function (j) { return px[j]; }), y: ai.map(function (j) { return W.ratio[j]; }), text: ai.map(tip), hoverinfo: 'text', marker: { color: sem.hit, size: 11, symbol: 'diamond' } });
      var L = baseLayout(t, { ylab: 'Observed / expected', legend: true, xaxis: winAxis(W, t),
        shapes: [{ type: 'line', xref: 'paper', x0: 0, x1: 1, yref: 'y', y0: 1, y1: 1, line: { color: CHROME[t].axis, width: 1, dash: 'dash' } }] });
      return { data: data, layout: L };
    });
    // z-score
    plot(z, function (t) {
      var sem = SEM[t], rx = [], ry = [], rt = [];
      W.ref_names.forEach(function (nm, r) { W.ref_z[r].forEach(function (v, j) { rx.push(px[j]); ry.push(v); rt.push(nm + '<br>exon ' + (W.xlab[j] || '')); }); });
      var ni = range(px.length).filter(function (j) { return !aff[j]; }), ai = range(px.length).filter(function (j) { return aff[j]; });
      var data = [
        { type: 'scatter', mode: 'markers', name: 'Reference samples', x: rx, y: ry, text: rt, hoverinfo: 'text+y', marker: { color: sem.ref, size: 5, opacity: 0.7 } },
        { type: 'scatter', mode: 'markers', name: 'Test sample', x: ni.map(function (j) { return px[j]; }), y: ni.map(function (j) { return W.z[j]; }), text: ni.map(function (j) { return 'exon ' + W.xlab[j]; }), hoverinfo: 'text+y', marker: { color: sem.test, size: 8 } },
        { type: 'scatter', mode: 'markers', name: 'Affected exon(s)', x: ai.map(function (j) { return px[j]; }), y: ai.map(function (j) { return W.z[j]; }), text: ai.map(function (j) { return 'exon ' + W.xlab[j]; }), hoverinfo: 'text+y', marker: { color: sem.hit, size: 11, symbol: 'diamond' } }
      ];
      var L = baseLayout(t, { ylab: 'Z-score', legend: true, xaxis: winAxis(W, t),
        shapes: [{ type: 'line', xref: 'paper', x0: 0, x1: 1, yref: 'y', y0: 0, y1: 0, line: { color: CHROME[t].axis, width: 1, dash: 'dash' } }] });
      return { data: data, layout: L };
    });
  }

  // ---- Cohort
  tab('cohort', 'Cohort', function () { return DATA.pca || DATA.corr || (DATA.depth_hist && DATA.depth_hist.counts.length); }, function (root) {
    var p1 = plotBox('tall'), p2 = plotBox('tall'), p3 = plotBox('short');
    var cards = [];
    if (DATA.pca) cards.push(card('PCA of coverage profiles', DATA.pca.n_samples + ' samples · ' + fmtInt(DATA.pca.n_rows) + ' targets (log2 counts)', p1));
    if (DATA.corr) cards.push(card('Sample–sample correlation', 'Spearman, targets of all chromosomes in the BED', p2));
    if (cards.length) root.appendChild(h('div', { class: 'grid g2' }, cards));
    if (!DATA.pca) root.appendChild(h('div', { style: { marginTop: '12px' } }, card('PCA of coverage profiles', null, emptyNote('PCA needs at least 3 samples with variable coverage.'))));
    if (DATA.depth_hist && DATA.depth_hist.counts.length)
      root.appendChild(h('div', { style: { marginTop: '12px' } }, card('Target depth distribution', 'median depth of each target across samples; dashed line = min_cov', p3)));
    setTimeout(function () { if (DATA.pca) pcaPlot(p1); if (DATA.corr) corrPlot(p2); if (DATA.depth_hist) depthHist(p3); }, 0);
  });
  function pcaPlot(div) {
    var P = DATA.pca;
    plot(div, function (t) {
      var groups = uniq(P.group), cc = catColors(Math.max(groups.length, 1), t), s = sel();
      var data = groups.map(function (g, k) {
        var ii = range(P.names.length).filter(function (i) { return P.group[i] === g; });
        return { type: 'scatter', mode: 'markers+text', name: g, x: ii.map(function (i) { return P.pc1[i]; }), y: ii.map(function (i) { return P.pc2[i]; }),
          text: ii.map(function (i) { return P.names[i]; }), textposition: 'top right', textfont: { size: 10, color: CHROME[t].secondary },
          hovertemplate: '<b>%{text}</b><br>PC1 %{x:.2f}<br>PC2 %{y:.2f}<extra>' + g + '</extra>',
          marker: { color: cc[k], size: ii.map(function (i) { return s && P.names[i] === s ? 15 : 9; }), line: { width: ii.map(function (i) { return s && P.names[i] === s ? 3 : 0; }), color: CHROME[t].text } } };
      });
      return { data: data, layout: baseLayout(t, { xlab: 'PC1 (' + fmt(P.var1) + '%)', ylab: 'PC2 (' + fmt(P.var2) + '%)', legend: groups.length > 1, margin: { l: 56, r: 24, t: 28, b: 48 } }) };
    });
  }
  function corrPlot(div) {
    var Cm = DATA.corr;
    plot(div, function (t) {
      var lo = Math.min.apply(null, Cm.z.map(function (r) { return Math.min.apply(null, r.filter(function (v) { return v !== null; })); }));
      var L = baseLayout(t, { xaxis: { type: 'category', tickangle: -45, showgrid: false }, yaxis: { type: 'category', autorange: 'reversed', showgrid: false }, margin: { l: 90, r: 16, t: 16, b: 90 } });
      return { data: [{ type: 'heatmap', x: Cm.names, y: Cm.names, z: Cm.z, zmin: Math.min(0.9, isFinite(lo) ? lo : 0.9), zmax: 1,
        colorscale: SEQ.map(function (c, i) { return [i / (SEQ.length - 1), c]; }), colorbar: { thickness: 10, len: 0.6, title: { text: 'r' } },
        hovertemplate: '%{y} × %{x}<br>r = %{z:.3f}<extra></extra>' }], layout: L };
    });
  }
  function depthHist(div) {
    var Hh = DATA.depth_hist;
    plot(div, function (t) {
      var sem = SEM[t];
      var data = [{ type: 'bar', x: Hh.mid, y: Hh.counts, marker: { color: Hh.mid.map(function (m) { return TH.min_cov !== null && TH.min_cov !== undefined && Math.pow(10, m) < TH.min_cov ? sem.bad : sem.dup; }) },
        hovertemplate: 'median depth ≈ %{customdata:.0f}<br>%{y} targets<extra></extra>', customdata: Hh.mid.map(function (m) { return Math.pow(10, m); }) }];
      var thr = Math.log10(Math.max(TH.min_cov || 1, 1));
      var hi = Math.ceil(Math.max.apply(null, Hh.mid) + 0.05), tv = range(hi + 1);
      var L = baseLayout(t, { xlab: 'Median depth per target (log scale)', ylab: 'Targets',
        xaxis: { type: 'linear', tickmode: 'array', tickvals: tv, ticktext: tv.map(function (e) { return fmtInt(Math.pow(10, e)); }) },
        shapes: TH.min_cov ? [{ type: 'line', xref: 'x', yref: 'paper', x0: thr, x1: thr, y0: 0, y1: 1, line: { color: sem.warn, width: 2, dash: 'dash' } }] : [] });
      return { data: data, layout: L };
    });
  }

  // ---- Methods
  tab('methods', 'Methods', function () { return true; }, function (root) {
    var blocks = [];
    var kv = h('dl', { class: 'kv' });
    DATA.settings.forEach(function (p) { kv.appendChild(h('dt', { text: p[0] })); kv.appendChild(h('dd', { text: p[1] })); });
    blocks.push(card('Settings', 'from the config used for this run', kv));
    var sc = DATA.scoring;
    if (sc) {
      var tbl = { cols: ['Category', 'Correlation', 'References', 'Fold change'], n: sc.rows.length, data: {} };
      tbl.cols.forEach(function (c, j) { tbl.data[c] = sc.rows.map(function (r) { return r[j]; }); });
      var hh = h('div');
      hh.appendChild(h('ul', null, h('li', null, h('b', { text: 'HIGH' }), ' only if all HIGH criteria are met.'), h('li', null, h('b', { text: 'MEDIUM' }), ' if HIGH is not met but all MEDIUM criteria are.'), h('li', null, h('b', { text: 'LOW' }), ' otherwise.'),
        sc.forced_low && sc.forced_low.length ? h('li', null, h('b', { text: 'Always LOW' }), ' (homology issues): ' + sc.forced_low.join(', ')) : null,
        sc.terminal ? h('li', { text: sc.terminal }) : null));
      var t = h('table', { class: 'tbl' }, h('thead', null, h('tr', null, tbl.cols.map(function (c) { return h('th', { text: c }); }))),
        h('tbody', null, sc.rows.map(function (r) { return h('tr', null, r.map(function (v, j) { return h('td', null, j === 0 ? chip(v, v) : v); })); })));
      hh.appendChild(h('div', { class: 'tbl-wrap' }, t));
      blocks.push(card('Confidence score', null, hh));
    }
    blocks.push(card('Reading the plots', null, h('div', null,
      h('p', { text: 'Each CNV call shows the test sample against its selected reference samples over the affected exons ± 5 neighbouring targets. Coverage is shown on a log scale after scaling each reference to the median coverage of the test sample in the window.' }),
      h('p', { text: 'Observed / expected is the test count divided by the count expected from the reference aggregate; the grey ribbon is the 95% interval of the beta-binomial model (overdispersion phi fitted by ExomeDepth). The z-score uses the model\'s own standard deviation, so reference samples should scatter around 0.'}),
      h('p', { text: 'If more background exons than expected fall outside the interval, the call is flagged "check" (poorly calibrated reference set).' }))));
    if (DATA.sample_table && DATA.sample_table.n) { var h1 = h('div'); tableView(h1, DATA.sample_table, { filename: 'echo_sample_table', pageSize: 15 }); blocks.push(card('Sample table', null, h1)); }
    if (DATA.ref_bams && DATA.ref_bams.n) { var h2 = h('div'); tableView(h2, DATA.ref_bams, { filename: 'echo_reference_bams', pageSize: 15 }); blocks.push(card('External reference BAMs', null, h2)); }
    root.appendChild(h('div', { class: 'stack' }, blocks));
  });

  // ---- Run
  tab('run', 'Run', function () { return true; }, function (root) {
    var L = DATA.log || {}, blocks = [];
    if (L.warnings && L.warnings.length) blocks.push(card('Warnings and errors from the log (' + L.warnings.length + ')', null, h('pre', { text: L.warnings.join('\n') })));
    if (DATA.config_yaml) blocks.push(card('Configuration used', null, h('pre', { text: DATA.config_yaml })));
    if (L.tail && L.tail.length) blocks.push(card('Pipeline log' + (L.truncated ? ' (tail)' : ''), null, h('pre', { text: L.tail.join('\n') })));
    if (DATA.session && DATA.session.length) blocks.push(card('Session info', null, h('pre', { text: DATA.session.join('\n') })));
    if (!blocks.length) blocks.push(emptyNote('No run information was recorded.'));
    root.appendChild(h('div', { class: 'stack' }, blocks));
  });

  // ------------------------------------------------------------------ shell
  var panels = {}, tabButtons = {}, active = null;
  function buildTab(t) {
    t.built = true; clear(panels[t.id]);
    try { t.build(panels[t.id]); } catch (e) { panels[t.id].appendChild(h('div', { class: 'empty', text: 'This tab could not be built: ' + e.message })); if (window.console) console.error(e); }
  }
  function activate(id, push) {
    var t = TABS.filter(function (x) { return x.id === id && x.available(); })[0] || TABS.filter(function (x) { return x.available(); })[0];
    if (!t) return;
    id = t.id;
    active = id;
    Object.keys(panels).forEach(function (k) {
      var on = k === id; panels[k].hidden = !on;
      tabButtons[k].setAttribute('aria-selected', on ? 'true' : 'false'); tabButtons[k].tabIndex = on ? 0 : -1;
    });
    if (!t.built) buildTab(t); else refreshVisible(panels[id]);
    setTimeout(function () { refreshVisible(panels[id]); }, 60);
    if (push) writeHash();
  }
  function writeHash() {
    try { history.replaceState(null, '', '#tab=' + active + (sel() && !META.focus ? '&sample=' + encodeURIComponent(sel()) : '')); } catch (e) { /* file:// may refuse */ }
  }
  function sampleChanged() {
    TABS.forEach(function (t) { if (t.built) { t.built = false; clear(panels[t.id]); } });
    registry = registry.filter(function (d) { return document.body.contains(d); });
    if (active) { var cur = TABS.filter(function (x) { return x.id === active; })[0]; buildTab(cur); setTimeout(function () { refreshVisible(panels[active]); }, 60); }
    writeHash(); updateBadges();
  }
  function updateBadges() {
    var box = $('#badges'); if (!box) return; clear(box);
    var nErr = activeAlerts(null, 'error').length, nWarn = activeAlerts(null, 'warning').length - nErr;
    box.appendChild(h('span', { class: 'badge', text: META.prefix ? META.prefix : 'ECHO' }));
    box.appendChild(h('span', { class: 'badge', text: 'modechrom ' + META.modechrom + ' · ' + META.modechrom_label }));
    if (META.genome_version) box.appendChild(h('span', { class: 'badge', text: META.genome_version }));
    if (META.focus) box.appendChild(h('span', { class: 'badge', text: 'single-sample report' }));
    if (nErr) box.appendChild(h('button', { class: 'badge bad', type: 'button', text: '✖ ' + plural(nErr, 'error'), title: 'Open QC & warnings', onclick: function () { nav('qc'); } }));
    if (nWarn > 0) box.appendChild(h('button', { class: 'badge warn', type: 'button', text: '▲ ' + plural(nWarn, 'warning'), title: 'Open QC & warnings', onclick: function () { nav('qc'); } }));
    if (!nErr && nWarn <= 0) box.appendChild(h('span', { class: 'badge ok', text: '✔ no warnings' }));
  }
  function init() {
    applyThemeAttr();
    document.title = META.title;
    $('#title').textContent = META.title;
    $('#generated').textContent = META.generated + (META.version ? ' · ECHO ' + META.version : '');
    var tabsEl = $('#tabs'), main = $('#panels');
    var avail = TABS.filter(function (t) { return t.available(); });
    avail.forEach(function (t, i) {
      var b = h('button', { class: 'tab', role: 'tab', id: 'tab-' + t.id, 'aria-controls': 'panel-' + t.id, 'aria-selected': 'false', tabindex: '-1', type: 'button', text: t.label,
        onclick: function () { activate(t.id, true); },
        onkeydown: function (e) {
          var k = e.key === 'ArrowRight' ? 1 : e.key === 'ArrowLeft' ? -1 : 0;
          if (!k) return; e.preventDefault();
          var nx = avail[(i + k + avail.length) % avail.length]; activate(nx.id, true); tabButtons[nx.id].focus();
        } });
      tabButtons[t.id] = b; tabsEl.appendChild(b);
      panels[t.id] = h('section', { class: 'panel', role: 'tabpanel', id: 'panel-' + t.id, 'aria-labelledby': 'tab-' + t.id, hidden: true });
      main.appendChild(panels[t.id]);
    });
    // sample selector (multi-sample reports only)
    var pick = $('#sample-pick');
    if (!META.focus && DATA.samples.n > 1) {
      var selEl = h('select', { id: 'sample-sel', 'aria-label': 'Sample' }, h('option', { value: '', text: 'All samples (' + DATA.samples.n + ')' }),
        DATA.samples.data.Sample.map(function (nm) { return h('option', { value: nm, text: nm }); }));
      selEl.addEventListener('change', function () { setSample(selEl.value); });
      pick.appendChild(h('label', null, 'Sample ')); pick.appendChild(selEl);
      pick.firstChild.setAttribute('for', 'sample-sel');
    } else if (pick) { pick.hidden = true; }
    onSample(sampleChanged);
    var btn = $('#theme-btn');
    var label = function () { btn.textContent = 'Theme: ' + themePref; };
    label();
    btn.addEventListener('click', function () {
      themePref = themePref === 'auto' ? 'light' : themePref === 'light' ? 'dark' : 'auto';
      store('echo-theme', themePref); applyThemeAttr(); label(); rerenderAll();
    });
    if (window.matchMedia) {
      var mq = window.matchMedia('(prefers-color-scheme: dark)');
      var onChange = function () { if (themePref === 'auto') rerenderAll(); };
      if (mq.addEventListener) mq.addEventListener('change', onChange);
    }
    window.addEventListener('beforeprint', function () {
      TABS.forEach(function (t) { if (t.available() && !t.built) buildTab(t); });
      refreshAll();
    });
    var mt = /tab=([\w-]+)/.exec(location.hash), ms = /sample=([^&]+)/.exec(location.hash);
    if (ms && !META.focus) { var want = decodeURIComponent(ms[1]); if (DATA.samples.data.Sample.indexOf(want) >= 0) setSample(want, true); }
    updateBadges();
    activate(mt ? mt[1] : 'overview', false);
    window.addEventListener('hashchange', function () { var mm = /tab=([\w-]+)/.exec(location.hash); if (mm && mm[1] !== active) activate(mm[1], false); });
    document.documentElement.setAttribute('data-ready', '1');
  }
  function refreshAll() { $$('.panel').forEach(function (p) { refreshVisible(p); }); }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
})();
