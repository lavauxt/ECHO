# ============================================================================
# Unified report assets: CSS + JavaScript + HTML assembly.
#
# Shared by the ECHO and CANOPE report generators so both packages emit
# byte-identical HTML structure, styling and client-side behaviour.  The only
# thing that differs between an ECHO and a CANOPE report is the JSON payload
# (the data blob embedded in <script id="unified-data">).
#
# Both echo_report.R and canope_report.R source() this file (or it is loaded
# once per package) and call report_html(payload, settings, out_file).
#
# Requires R >= 4.0 for raw strings (r"---(...)---").
# ============================================================================

# ----------------------------------------------------------------------------
# Stylesheet  (merged ECHO + CANOPE design system)
# ----------------------------------------------------------------------------

REPORT_CSS <- r"---(/* Unified interactive report -- design tokens, layout and components.
   Merges the ECHO and CANOPE stylesheets so both packages share one look.
   Categorical hues are applied by JS in fixed order; chrome (surface / ink /
   grid) swaps between light and dark.  Colour is never the only cue: every
   status colour is paired with an icon, symbol or text label. */
:root {
  color-scheme: light;
  --page: #f9f9f7;
  --surface: #fcfcfb;
  --text: #0b0b0b;
  --text-2: #52514e;
  --muted: #898781;
  --grid: #e1e0d9;
  --axis: #c3c2b7;
  --ring: rgba(11, 11, 11, 0.10);
  --accent: #2a78d6;
  --accent-ink: #ffffff;
  --accent-soft: #e7f0fc;
  --good: #006300;
  --warn: #8a5a00;
  --fail: #b3261e;
  --bad: #c42b2b;
  --ok-ink: #0b6b3a;
  --fail-soft: #fdeceb;
  --warn-soft: #fff3d6;
  --ok-soft: #e4f5ea;
  --info-soft: #e7f0fc;
  --bad-soft: #fdeaea;
  --c-del: #e34948; --c-dup: #2a78d6;
  --c-high: #1baf7a; --c-medium: #eda100; --c-low: #898781; --c-na: #c3c2b7;
  --radius: 10px;
  --font: system-ui, -apple-system, "Segoe UI", sans-serif;
  --mono: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
}
@media (prefers-color-scheme: dark) {
  :root:where(:not([data-theme="light"])) {
    color-scheme: dark;
    --page: #0d0d0d; --surface: #1a1a19; --text: #f0efec; --text-2: #c3c2b7;
    --muted: #898781; --grid: #2c2c2a; --axis: #383835;
    --ring: rgba(255, 255, 255, 0.10); --accent: #3987e5; --accent-ink: #0b0b0b;
    --accent-soft: #1c2a3d; --good: #0ca30c; --warn: #fab219;
    --fail: #ff8a80; --bad: #f08080; --ok-ink: #5fd38d;
    --fail-soft: #3a1c1a; --warn-soft: #3a2d0b; --ok-soft: #12301f;
    --info-soft: #1c2a3d; --bad-soft: #3a1c1c;
    --c-del: #e66767; --c-dup: #3987e5; --c-high: #199e70; --c-medium: #c98500; --c-low: #898781; --c-na: #555552;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --page: #0d0d0d; --surface: #1a1a19; --text: #f0efec; --text-2: #c3c2b7;
  --muted: #898781; --grid: #2c2c2a; --axis: #383835;
  --ring: rgba(255, 255, 255, 0.10); --accent: #3987e5; --accent-ink: #0b0b0b;
  --accent-soft: #1c2a3d; --good: #0ca30c; --warn: #fab219;
  --fail: #ff8a80; --bad: #f08080; --ok-ink: #5fd38d;
  --fail-soft: #3a1c1a; --warn-soft: #3a2d0b; --ok-soft: #12301f;
  --info-soft: #1c2a3d; --bad-soft: #3a1c1c;
  --c-del: #e66767; --c-dup: #3987e5; --c-high: #199e70; --c-medium: #c98500; --c-low: #898781; --c-na: #555552;
}

* { box-sizing: border-box; }
html { scroll-behavior: smooth; }
body {
  margin: 0; background: var(--page); color: var(--text);
  font: 14px/1.5 var(--font); -webkit-font-smoothing: antialiased;
}
a { color: var(--accent); }
h1, h2, h3 { margin: 0; font-weight: 600; letter-spacing: -0.01em; }
h2 { font-size: 18px; } h3 { font-size: 14px; }
code, pre { font-family: var(--mono); font-size: 12px; }
pre { margin: 0; padding: 12px; overflow: auto; background: var(--page);
      border: 1px solid var(--ring); border-radius: 8px; max-height: 420px; }
.muted { color: var(--muted); } .small { font-size: 12px; }
.sr-only { position: absolute; width: 1px; height: 1px; overflow: hidden; clip: rect(0 0 0 0); }

/* ---- Header & tabs ---- */
.top { position:sticky; top:0; z-index:100; background:var(--surface); border-bottom:1px solid var(--ring); box-shadow:0 2px 10px var(--ring); }
.top-inner { max-width: 1440px; margin: 0 auto; padding: 14px 16px 0; }
.title-row { display:flex; align-items:center; column-gap:24px; row-gap:14px; flex-wrap:wrap; }
.title-row h1 { font-size: 20px; }
.badge { display: inline-flex; align-items: center; gap: 6px; padding: 2px 10px;
         border-radius: 999px; background: var(--accent-soft); color: var(--text-2);
         font-size: 12px; border: 1px solid var(--ring); }
.badge.demo { background: transparent; border-color: var(--warn); color: var(--warn); font-weight: 600; }
.badge.warn { background: var(--warn-soft); border-color: var(--warn); color: var(--warn); font-weight: 600; cursor: pointer; }
.badge.bad, .badge.fail { background: var(--fail-soft); border-color: var(--fail); color: var(--fail); font-weight: 600; cursor: pointer; }
.badge.ok { color: var(--good); }
.spacer { flex: 1; }
.tabs { display: flex; gap: 2px; overflow-x: auto; scrollbar-width: thin; margin-top: 10px; }
.tab { appearance: none; border: 0; background: none; color: var(--text-2); font: inherit;
       padding: 10px 14px; cursor: pointer; border-bottom: 2px solid transparent;
       white-space: nowrap; border-radius: 6px 6px 0 0; }
.tab:hover { background: var(--accent-soft); color: var(--text); }
.tab[aria-selected="true"] { color: var(--text); border-bottom-color: var(--accent); font-weight: 600; }
.tab:focus-visible, button:focus-visible, select:focus-visible, input:focus-visible {
  outline: 2px solid var(--accent); outline-offset: 2px; }
main { max-width: 1440px; margin: 0 auto; padding: 16px; }
.panel[hidden] { display: none; }

/* ---- Cards, grids, controls ---- */
.collapsible-section { padding: 0; overflow: hidden; }
.collapsible-section > summary { cursor:pointer; padding:14px 16px; font-weight:650; list-style:none; }
.collapsible-section > summary::-webkit-details-marker { display:none; }
.collapsible-section > summary::after { content:'＋'; float:right; color:var(--text-2); }
.collapsible-section[open] > summary::after { content:'−'; }
.collapsible-section > .stack { padding:0 14px 14px; }
.card { background: var(--surface); border: 1px solid var(--ring); border-radius: var(--radius);
        padding: 14px; min-width: 0; }
.card > header { display: flex; align-items: baseline; gap: 8px; margin-bottom: 8px; flex-wrap: wrap; }
.card > header .hint { margin-left: auto; }
.card.flagged { border-color: var(--warn); }
.grid { display: grid; gap: 12px; }
.g2 { grid-template-columns: repeat(auto-fit, minmax(min(100%, 460px), 1fr)); }
.g3 { grid-template-columns: repeat(auto-fit, minmax(min(100%, 320px), 1fr)); }
.g-main { grid-template-columns: minmax(0, 1fr) 360px; align-items: start; }
.g-viewer { grid-template-columns: 340px minmax(0, 1fr); align-items: start; }
.stack { display: flex; flex-direction: column; gap: 12px; }
.stack.tight { gap: 6px; }
.toolbar { display: flex; gap: 12px; flex-wrap: wrap; align-items: flex-end; margin-bottom: 12px; }
.toolbar:last-child { margin-bottom: 0; }
.field { display: flex; flex-direction: column; gap: 3px; font-size: 12px; color: var(--text-2); min-width: 0; max-width: 100%; }
.field.inline { flex-direction: row; align-items: center; gap: 6px; }
.field.inline input[type="checkbox"] { width: 16px; height: 16px; accent-color: var(--accent); }
select, input[type="text"], input[type="number"], input[type="search"] {
  font: inherit; color: var(--text); background: var(--surface); border: 1px solid var(--axis);
  border-radius: 6px; padding: 5px 8px; min-width: 120px; max-width: 280px; }
input[type="range"] { accent-color: var(--accent); }
button.btn { font: inherit; padding: 5px 12px; border-radius: 6px; cursor: pointer;
             border: 1px solid var(--axis); background: var(--surface); color: var(--text); }
button.btn:hover { background: var(--accent-soft); }
button.btn.primary { background: var(--accent); border-color: var(--accent); color: var(--accent-ink); }
button.btn[disabled] { opacity: 0.45; cursor: not-allowed; }
button.btn.small-btn { padding: 1px 9px; font-size: 12px; margin-left: 6px; }
.plot { width: 100%; min-height: 360px; }
.plot.tall { min-height: 560px; } .plot.short { min-height: 280px; }
.plot.strip { min-height: 70px; height: 70px; }
.plot.mid { min-height: 240px; }
.plot.viewer { height: 760px; }
.empty { padding: 28px 12px; text-align: center; color: var(--muted); }

/* ---- Stat tiles ---- */
.tiles { display: grid; gap: 12px; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); }
.tile { background: var(--surface); border: 1px solid var(--ring); border-radius: var(--radius); padding: 12px 14px; }
.tile .v { font-size: 26px; font-weight: 600; letter-spacing: -0.02em; }
.tile .l { color: var(--text-2); font-size: 12px; }
.tile .detail { display:block; margin-top:3px; color:var(--muted); font-size:11px; }
.tile.ok { border: 1px solid var(--c-high); border-left: 4px solid var(--c-high); }
.tile.warn { border: 1px solid var(--warn); border-left: 4px solid var(--warn); background: var(--warn-soft); }
.tile.fail, .tile.bad { border: 1px solid var(--fail); border-left: 4px solid var(--fail); background: var(--fail-soft); }
.tile.warn .v { color: var(--warn); } .tile.fail .v, .tile.bad .v { color: var(--fail); } .tile.ok .v { color: var(--ok-ink); }
.scope { display: grid; gap: 12px; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); margin-bottom: 12px; }
.scope .tile .l { margin-top: 2px; }

/* ---- Legend (HTML so categories can be toggled) ---- */
.legend { display: flex; flex-wrap: wrap; gap: 4px 12px; margin-top: 8px; font-size: 12px; color: var(--text-2); }
.legend button { appearance: none; border: 0; background: none; font: inherit; color: inherit;
                 display: inline-flex; gap: 6px; align-items: center; cursor: pointer; padding: 2px 4px; border-radius: 4px; }
.legend button:hover { background: var(--accent-soft); }
.legend button[aria-pressed="false"] { opacity: 0.4; text-decoration: line-through; }
.swatch { width: 10px; height: 10px; border-radius: 3px; display: inline-block; }
.swatch.conf.HIGH { background: var(--c-high); } .swatch.conf.MEDIUM { background: var(--c-medium); }
.swatch.conf.LOW { background: var(--c-low); } .swatch.conf.NA { background: var(--c-na); }

/* ---- Tables ---- */
.tbl-wrap { overflow: auto; border: 1px solid var(--ring); border-radius: 8px; max-height: 520px; }
table.tbl { border-collapse: collapse; width: 100%; font-size: 12.5px; }
.tbl th, .tbl td { padding: 5px 10px; text-align: left; white-space: nowrap; border-bottom: 1px solid var(--grid); }
.tbl th { position: sticky; top: 0; background: var(--surface); cursor: pointer; user-select: none;
          font-weight: 600; color: var(--text-2); border-bottom: 1px solid var(--axis); z-index: 1; }
.tbl th[aria-sort="ascending"]::after { content: " \25B2"; font-size: 9px; }
.tbl th[aria-sort="descending"]::after { content: " \25BC"; font-size: 9px; }
.tbl td.num { text-align: right; font-variant-numeric: tabular-nums; }
.tbl td.wrap { white-space: normal; min-width: 260px; max-width: 520px; }
.tbl tbody tr:hover { background: var(--accent-soft); }
.tbl tbody tr.clickable, .tbl tbody tr.pick { cursor: pointer; }
.tbl tbody tr.sel { outline: 2px solid var(--accent); outline-offset: -2px; background: var(--accent-soft); }
.tbl td.warn { background: var(--warn-soft); color: var(--warn); font-weight: 600; }
.tbl td.fail { background: var(--fail-soft); color: var(--fail); font-weight: 600; }
.tbl-bar { display: flex; gap: 10px; align-items: center; margin: 8px 0; flex-wrap: wrap; }
.pager { display: inline-flex; gap: 6px; align-items: center; margin-left: auto; }

/* ---- Selection bars / kv / notes ---- */
.bars { display: flex; flex-direction: column; gap: 5px; }
.bar-row { display: grid; grid-template-columns: 110px 1fr 70px; gap: 8px; align-items: center; font-size: 12px; }
.bar-row .lab { overflow: hidden; text-overflow: ellipsis; white-space: nowrap; color: var(--text-2); }
.bar-track { background: var(--grid); border-radius: 4px; height: 10px; overflow: hidden; }
.bar-fill { height: 100%; border-radius: 4px; }
.bar-fill.lv-ok { background: var(--c-high); }
.bar-fill.lv-warn { background: var(--c-medium); }
.bar-fill.lv-fail { background: var(--c-del); }
.note { font-size: 12px; color: var(--muted); }
.note-flag { font-size: 11.5px; }
.kv { display: grid; grid-template-columns: max-content 1fr; gap: 3px 14px; font-size: 12.5px; }
.kv dt { color: var(--text-2); } .kv dd { margin: 0; overflow-wrap: anywhere; }
.prose { width:100%; max-width:none; }
.prose p { margin: 0 0 10px; width:100%; max-width:none; }

/* ---- Alerts (always an icon + text, never colour alone) ---- */
.alerts { display: flex; flex-direction: column; gap: 8px; margin-bottom: 12px; }
.alert { display: flex; gap: 10px; align-items: flex-start; padding: 9px 14px; border-radius: 8px;
         border: 1px solid var(--ring); border-left: 4px solid var(--muted); background: var(--surface); font-size: 13px; }
.alert .ico, .alert .ic { font-weight: 700; flex: none; width: 18px; text-align: center; line-height: 1.35; }
.alert .body { min-width: 0; overflow-wrap: anywhere; }
.alert .ttl { font-weight: 600; }
.alert.fail, .alert.error { border: 1px solid var(--fail); border-left: 4px solid var(--fail); background: var(--fail-soft); }
.alert.fail .ico, .alert.fail .ic, .alert.fail .ttl, .alert.error .ico, .alert.error .ic, .alert.error .ttl { color: var(--fail); }
.alert.warn, .alert.warning { border: 1px solid var(--warn); border-left: 4px solid var(--warn); background: var(--warn-soft); }
.alert.warn .ico, .alert.warn .ic, .alert.warn .ttl, .alert.warning .ico, .alert.warning .ic, .alert.warning .ttl { color: var(--warn); }
.alert.ok { border: 1px solid var(--c-high); border-left: 4px solid var(--c-high); background: var(--ok-soft); }
.alert.ok .ico, .alert.ok .ic { color: var(--ok-ink); }
.alert.info { border: 1px solid var(--accent); border-left: 4px solid var(--accent); background: var(--info-soft); }
.alert.info .ico, .alert.info .ic, .alert.info .ttl { color: var(--accent); }
.alert a { cursor: pointer; }

/* Alert content is controlled by the single global severity selector. */

/* ---- Chips ---- */
.chip { display: inline-flex; align-items: center; gap: 4px; padding: 0 8px; border-radius: 999px; font-size: 11.5px;
        font-weight: 600; line-height: 20px; border: 1px solid var(--ring); white-space: nowrap; background: var(--page); color: var(--text-2); }
.chip .ico { font-size: 11px; }
.chip.ok { color: var(--ok-ink); background: var(--ok-soft); border-color: var(--c-high); }
.chip.warn, .chip.warning { color: var(--warn); background: var(--warn-soft); border-color: var(--warn); }
.chip.fail, .chip.error { color: var(--fail); background: var(--fail-soft); border-color: var(--fail); }
.chip.info { color: var(--accent); background: var(--info-soft); border-color: var(--accent); }
.chip.cat { font-weight: 500; }
.chip.type.DEL, .chip.DEL { color: var(--c-del); border-color: var(--c-del); background: transparent; }
.chip.type.DUP, .chip.DUP { color: var(--c-dup); border-color: var(--c-dup); background: transparent; }
.chip.conf.HIGH, .chip.HIGH { color: var(--ok-ink); background: var(--ok-soft); border-color: var(--c-high); }
.chip.conf.MEDIUM, .chip.MEDIUM { color: var(--warn); background: var(--warn-soft); border-color: var(--warn); }
.chip.conf.LOW, .chip.LOW { color: var(--text-2); background: var(--page); border-color: var(--muted); }
.chip.conf.NA, .chip.NA { color: var(--text-2); background: var(--page); border-color: var(--axis); }
.chip-btn { appearance: none; font: inherit; font-size: 11.5px; font-weight: 600; display: inline-flex; gap: 4px; align-items: center;
            padding: 0 8px; line-height: 20px; border-radius: 999px; cursor: pointer; border: 1px solid var(--ring); background: var(--page); color: var(--text); }
.chip-btn:hover { background: var(--accent-soft); }
.chip-btn.warn { color: var(--warn); background: var(--warn-soft); }
.chip-btn.fail { color: var(--fail); background: var(--fail-soft); }
.chip-btn.link { color: var(--accent); }
.chips { display: flex; flex-wrap: wrap; gap: 6px; }
.count-badge { display: inline-block; margin-left: 6px; min-width: 18px; padding: 0 6px; border-radius: 999px; font-size: 11px; line-height: 18px;
               text-align: center; background: var(--warn); color: var(--page); font-weight: 700; }
.pill-list { display: flex; flex-wrap: wrap; gap: 6px; }

/* ---- Segmented control ---- */
.seg { display: inline-flex; border: 1px solid var(--axis); border-radius: 6px; overflow: hidden; }
.seg button { font: inherit; border: 0; background: var(--surface); color: var(--text-2); padding: 4px 10px; cursor: pointer;
              border-right: 1px solid var(--axis); }
.seg button:last-child { border-right: 0; }
.seg button[aria-pressed="true"] { background: var(--accent); color: var(--accent-ink); }

/* ---- Sub-tabs (inner switch within a tab) ---- */
.subtabs { display: flex; gap: 2px; border-bottom: 1px solid var(--ring); margin-bottom: 14px; overflow-x: auto; scrollbar-width: thin; }
.subtab { appearance: none; border: 0; background: none; color: var(--text-2); font: inherit; font-size: 13px;
          padding: 8px 14px; cursor: pointer; border-bottom: 2px solid transparent; white-space: nowrap; border-radius: 6px 6px 0 0; }
.subtab:hover { background: var(--accent-soft); color: var(--text); }
.subtab[aria-selected="true"] { color: var(--text); border-bottom-color: var(--accent); font-weight: 600; }
.subpanel[hidden] { display: none; }

/* ---- Header filters ---- */
.sample-pick { display:flex; align-items:center; column-gap:22px; row-gap:12px; flex-wrap:wrap; font-size:12px; color:var(--text-2); }
.sample-pick select { width:min(260px, 30vw); min-width:170px; }
.sample-pick > label { font-weight:650; white-space:nowrap; }
.global-alert-filter { display:inline-flex; align-items:center; gap:10px; margin-left:4px; padding-left:20px; border-left:1px solid var(--ring); white-space:nowrap; }
.global-notice { max-width:100%; padding:0 0 8px; }
.global-notice:empty { display:none; }
.global-notice .alert { margin:0; padding:7px 12px; font-size:12.5px; }
.global-notice .alert .btn { margin-left:8px; padding:2px 9px; font-size:11.5px; }
.global-alert-filter > span { font-weight:650; font-size:12px; color:var(--text-2); }
.global-alert-filter .seg button { padding:5px 10px; font-size:11.5px; }
#theme-btn { margin-left:10px; }
.alert[hidden] { display:none !important; }

/* ---- Independently collapsible Methods and Run blocks ---- */
.collapsible-section { padding:0; overflow:hidden; }
.collapsible-section > summary { cursor:pointer; padding:12px 15px; font-weight:650; list-style:none; display:flex; align-items:baseline; gap:10px; }
.collapsible-section > summary::-webkit-details-marker { display:none; }
.collapsible-section > summary::after { content:'＋'; margin-left:auto; color:var(--text-2); }
.collapsible-section[open] > summary::after { content:'−'; }
.collapsible-section .collapse-hint { color:var(--muted); font-size:11.5px; font-weight:400; }
.collapsible-body { padding:0 14px 14px; }
@media (max-width:760px) {
  .sample-pick { column-gap:12px; row-gap:9px; }
  .sample-pick select { width:min(220px, 65vw); min-width:145px; }
  .global-alert-filter { margin-left:0; padding-left:12px; gap:6px; flex-wrap:wrap; }
  .global-alert-filter .seg button { padding:4px 6px; }
  #theme-btn { margin-left:0; }
}

/* ---- Chromosome strip ---- */
.chrom-strip { display: grid; gap: 10px; grid-template-columns: repeat(auto-fit, minmax(210px, 1fr)); }
.chrom-box { border: 1px solid var(--ring); border-left: 4px solid var(--muted); border-radius: 8px; padding: 8px 12px; display: flex; flex-direction: column; gap: 4px; align-items: flex-start; }
.chrom-box .nm { font-weight: 600; }
.chrom-box.ok { border-color: var(--c-high); border-left: 4px solid var(--c-high); }
.chrom-box.warn { border-color: var(--warn); border-left: 4px solid var(--warn); background: var(--warn-soft); }
.chrom-box.fail { border-color: var(--fail); border-left: 4px solid var(--fail); background: var(--fail-soft); }
.chrom-box.info { border-color: var(--accent); border-left: 4px solid var(--accent); }

/* ---- Gallery (per-call plots) ---- */
.gallery { display: grid; gap: 12px; grid-template-columns: repeat(auto-fill, minmax(280px, 1fr)); }
.shot { background: var(--surface); border: 1px solid var(--ring); border-radius: 8px; overflow: hidden; margin: 0; }
.shot button { all: unset; display: block; cursor: zoom-in; width: 100%; }
.shot img { display: block; width: 100%; height: 190px; object-fit: contain; background: #fff; }
.shot figcaption { padding: 6px 10px; font-size: 12px; color: var(--text-2); overflow-wrap: anywhere; }
.shot .nopreview { height: 190px; display: grid; place-items: center; color: var(--muted); background: var(--page); }
dialog.lightbox { border: 0; padding: 0; background: transparent; max-width: 96vw; max-height: 96vh; }
dialog.lightbox::backdrop { background: rgba(0, 0, 0, 0.7); }
dialog.lightbox img { max-width: 96vw; max-height: 92vh; background: #fff; border-radius: 8px; display: block; }

/* ---- Call viewer detail ---- */
.detail-head { display: flex; flex-wrap: wrap; gap: 8px 14px; align-items: baseline; margin-bottom: 8px; }
.detail-head h2 { font-size: 16px; }

@media (max-width: 1000px) { .g-main, .g-viewer, .g2, .g3 { grid-template-columns: minmax(0, 1fr); } }
@media (max-width: 720px) { main { padding: 12px; } .bar-row { grid-template-columns: 80px 1fr 56px; } }


/* ---- Overview: unified cohort summary ---- */
.overview-summary { display:grid; grid-template-columns:minmax(0, 1fr) minmax(0, 1.15fr); gap:18px; align-items:stretch; padding:16px; border:1px solid var(--ring); border-radius:var(--radius); background:var(--page); }
.overview-summary-single { grid-template-columns:minmax(0, 1fr); }
.overview-summary-single .overview-summary-stats, .overview-summary-single .overview-not-analysed { grid-column:1 / -1; }
.overview-summary-status, .overview-summary-stats { min-width:0; display:flex; flex-direction:column; }
.overview-summary-status { padding-right:16px; border-right:1px solid var(--ring); }
.overview-section-label { font-size:12px; font-weight:700; color:var(--text-2); margin:0 0 12px; }
.overview-summary .health-tile-grid { grid-template-columns:repeat(3, minmax(0, 1fr)); align-content:center; flex:1; }
.overview-summary .health-tile { min-width:0; min-height:112px; }
.overview-summary .tiles { grid-template-columns:repeat(3, minmax(0, 1fr)); grid-auto-rows:1fr; align-content:stretch; flex:1; }
.overview-summary .tile { min-width:0; padding:10px; }
.overview-summary .tile .v { font-size:22px; overflow-wrap:anywhere; }
.overview-summary .tile .l { font-size:11px; line-height:1.3; }
.overview-not-analysed { margin-top:12px; padding:11px 12px; border:1px solid var(--fail); border-left:4px solid var(--fail); border-radius:8px; background:var(--fail-bg); color:var(--text); }
.overview-not-analysed strong { display:block; margin-bottom:5px; color:var(--fail); font-size:12px; }
.overview-not-analysed ul { margin:0; padding-left:18px; display:flex; flex-wrap:wrap; gap:3px 22px; font-size:12px; }
.global-notice { max-width:100%; padding:8px 0 14px; }
.global-notice .alert { margin:0; padding:12px 16px; min-height:58px; }
.global-notice .alert .btn { margin-left:14px; padding:4px 11px; }
@media (max-width:900px) { .overview-summary { grid-template-columns:minmax(0, 1fr); } .overview-summary-status { padding:0 0 16px; border-right:0; border-bottom:1px solid var(--ring); } }
@media (max-width:560px) { .overview-summary .tiles { grid-template-columns:repeat(2, minmax(0, 1fr)); } }
/* Gene labels use an HTML track so long symbols remain readable at every width. */
.gene-track-container { width:100%; min-width:0; padding:10px 12px; overflow-x:auto; border:1px solid var(--ring); border-radius:8px; background:var(--surface); }
.gene-track { display:flex; align-items:stretch; gap:7px; width:100%; min-width:max-content; }
.gene-track-segment { flex:1 1 0; min-width:125px; min-height:62px; display:flex; flex-direction:column; justify-content:center; align-items:center; gap:4px; padding:8px 10px; border:1px solid var(--track-border, var(--accent)); border-bottom:3px solid var(--track-border, var(--accent)); border-radius:7px; text-align:center; background:var(--track-bg, var(--accent-soft)); }
.gene-track-segment.off { border-style:dashed; border-color:var(--axis); background:var(--page); }
.gene-track-name { font-size:12.5px; font-weight:700; line-height:1.25; overflow-wrap:anywhere; color:var(--text); }
.gene-track-meta { font-size:11px; color:var(--text-2); }
.gene-track-segment.track-a { --track-border:var(--accent); --track-bg:var(--accent-soft); }
.gene-track-segment.track-b { --track-border:var(--c-high); --track-bg:var(--ok-soft); }
.gene-track-segment.track-c { --track-border:var(--c-medium); --track-bg:var(--warn-soft); }
.gene-track-segment.track-d { --track-border:var(--c-del); --track-bg:var(--bad-soft); }
.gene-track-segment.track-e { --track-border:var(--c-dup); --track-bg:var(--info-soft); }
.gene-track-segment.track-f { --track-border:var(--muted); --track-bg:var(--page); }

/* ---- Overview: compact square QC status tiles ---- */
.health-tile-grid { display:grid; grid-template-columns:repeat(3, minmax(0, 1fr)); gap:12px; }
.health-tile { min-width:0; min-height:78px; display:flex; flex-direction:column; justify-content:center; align-items:center; text-align:center; padding:10px; border-radius:var(--radius); }
.health-tile .v { font-size:22px; font-weight:600; }
.health-tile .l { font-size:11px; line-height:1.3; font-weight:650; }
.health-tile .detail { font-size:11px; color:var(--muted); margin-top:3px; }
.calls-toolbar { align-items:flex-end; }
.calls-toolbar .field { min-width:max-content; }
.calls-toolbar .seg button { padding:4px 8px; font-size:11px; }
@media (max-width:560px) {
  .health-tile-grid { grid-template-columns:repeat(3, minmax(0, 1fr)); }
  .health-tile { min-height:98px; }
}

@media print {
  .top { position:static; box-shadow:none; }
  .tabs, .toolbar, .btn, .subtabs, .sample-pick, .seg { display: none !important; }
  .panel[hidden], .subpanel[hidden] { display: block !important; }
  .collapsible-section:not([open]) > :not(summary) { display: block !important; }
  .collapsible-section > summary { cursor: default; padding: 8px 0; }
  .collapsible-section > summary::after { content: ''; }
  .card { break-inside: avoid; }
}
@media (forced-colors: active) { .bar-fill, .swatch { forced-color-adjust: none; } }
)---"


# ----------------------------------------------------------------------------
# Client-side application  (unified -- handles both ECHO and CANOPE payloads)
# ----------------------------------------------------------------------------

REPORT_JS <- r"---(/* Unified interactive report -- client-side application.
   Reads the JSON blob in #unified-data and builds lazily-rendered tabs:
   Overview (Samples / Chromosomes / Call landscape) · Samples · QC & warnings · CNV calls · Methods · Run. The header contains one sample filter and one global severity filter. The payload's meta.source field ("ECHO" or "CANOPE") selects
   which source-specific call columns and panels are shown.  The ECHO call table
   follows the legacy report's explicit fields; CANOPE exposes its own model fields.

   Architecture follows the SPARKS / ECHO / CANOPE report family: h() DOM
   helper, tableView, theme switch, hash navigation, independently collapsible Methods and Run cards, a compact global severity filter,
   direct QC subtabs and a global sample filter.  No server and no network access are
   needed once plotly.js is loaded. */
(function () {
  'use strict';

  // ================================================================ data & utils
  var DATA = JSON.parse(document.getElementById('unified-data').textContent);
  var META = DATA.meta || {};
  var IS_ECHO = META.source === 'ECHO';
  var IS_CANOPE = META.source === 'CANOPE';
  var TH = DATA.thresholds || {};
  var SUM = DATA.summary || {};
  var S = DATA.samples || { cols: [], data: {}, n: 0 };
  var CALLS = DATA.calls || { cols: [], data: {}, n: 0 };
  var PLOTS = DATA.plots || DATA.windows || {};
  var CHR = DATA.chromosomes || {};
  var SCOPE = DATA.scope || { cols: [], data: {}, n: 0 };
  var QC = DATA.qc || null;             // CANOPE QC metrics block
  var QC_FLAGS = DATA.qc_flags || null;  // ECHO QC flags table
  var PCA = DATA.pca || null;
  var DEPTH_HIST = DATA.depth_hist || null;
  var CONF_RULES = DATA.conf_rules || {};
  var PARAMS = DATA.params || {};
  var SCORING = DATA.scoring || null;
  var SETTINGS = DATA.settings || [];
  var LOG = DATA.log || {};

  var SAMPLES = (S.data.sample || S.data.Sample || []).slice();
  var SAMPLE_KEY = S.data.sample ? 'sample' : 'Sample';
  var SIDX = {}; SAMPLES.forEach(function (s, i) { SIDX[s] = i; });

  var $ = function (s, r) { return (r || document).querySelector(s); };
  var $$ = function (s, r) { return Array.prototype.slice.call((r || document).querySelectorAll(s)); };

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
  function uniq(a) { var seen = {}; return a.filter(function (v) { var k = String(v); if (seen[k]) return false; seen[k] = true; return true; }); }
  function asArray(x) { return x === null || x === undefined ? [] : Array.isArray(x) ? x : [x]; }
  function isNum(v) { return typeof v === 'number' && isFinite(v); }
  function fmtInt(n) { return Number(n).toLocaleString('en-US'); }
  function fmt(v) {
    if (v === null || v === undefined || (typeof v === 'number' && !isFinite(v))) return '';
    if (typeof v === 'boolean') return v ? 'yes' : 'no';
    if (typeof v !== 'number') return String(v);
    if (Number.isInteger(v)) return String(v);
    if (v !== 0 && Math.abs(v) < 1e-3) return v.toExponential(2);
    return String(Number(v.toPrecision(4)));
  }
  function orNa(v) { return (v === null || v === undefined || v === '' || (typeof v === 'number' && !isFinite(v))) ? 'n/a' : fmt(v); }
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
    try { if (val === undefined) return window.localStorage.getItem(key); window.localStorage.setItem(key, val); } catch (e) {}
    return null;
  }
  function plural(n, w, ws) { return n + ' ' + (n === 1 ? w : (ws || w + 's')); }
  function pretty(c) { return String(c).replace(/_/g, ' '); }

  // table helper: get a column array from a table object
  function tcol(T, c) { return (T && T.data && T.data[c]) || []; }
  function sval(s, col) { var i = SIDX[s]; return i === undefined ? null : tcol(S, col)[i]; }
  function hasCol(T, c) { return T && T.data && T.data[c] && T.data[c].length; }

  // ---- icons & levels ----
  var ICON = { fail: '\u2716', error: '\u2716', warn: '\u25b2', warning: '\u25b2', info: '\u2139', ok: '\u2714' };
  var LEVEL_LABEL = { fail: 'Fail', error: 'Fail', warn: 'Warning', warning: 'Warning', info: 'Note', ok: 'OK' };

  // ================================================================ theme & palette
  var CHROME = {
    light: { surface: '#fcfcfb', text: '#0b0b0b', secondary: '#52514e', muted: '#898781', grid: '#e1e0d9', axis: '#c3c2b7', mid: '#f0efec', band: '#f0efec' },
    dark: { surface: '#1a1a19', text: '#f0efec', secondary: '#c3c2b7', muted: '#898781', grid: '#2c2c2a', axis: '#383835', mid: '#383835', band: '#232321' }
  };
  var CAT = {
    light: ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948'],
    dark: ['#3987e5', '#d95926', '#199e70', '#c98500', '#d55181', '#008300', '#9085e9', '#e66767']
  };
  var SEM = {
    light: { DEL: '#e34948', DUP: '#2a78d6', HIGH: '#006300', MEDIUM: '#eda100', LOW: '#898781', NA: '#c3c2b7',
             fail: '#b3261e', error: '#b3261e', warn: '#8a5a00', warning: '#8a5a00', ok: '#1baf7a', info: '#2a78d6',
             test: '#2a78d6', ref: '#a9a8a0', flank: '#898781', hit: '#e34948', del: '#e34948', dup: '#2a78d6',
             high: '#006300', medium: '#eda100', low: '#898781', bad: '#c42b2b' },
    dark: { DEL: '#e66767', DUP: '#3987e5', HIGH: '#0ca30c', MEDIUM: '#fab219', LOW: '#898781', NA: '#555552',
            fail: '#ff8a80', error: '#ff8a80', warn: '#fab219', warning: '#fab219', ok: '#199e70', info: '#3987e5',
            test: '#3987e5', ref: '#6f6e68', flank: '#898781', hit: '#e66767', del: '#e66767', dup: '#3987e5',
            high: '#0ca30c', medium: '#fab219', low: '#898781', bad: '#f08080' }
  };
  var SEQ = ['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b'];
  var DIV = ['#e34948', '#f0a3a2', '#f0efec', '#9ec5f4', '#2a78d6'];
  var FONT = "system-ui, -apple-system, 'Segoe UI', sans-serif";

  var themePref = store('unified-report-theme') || 'auto';
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
  function alpha(hex, a) {
    var v = [1, 3, 5].map(function (i) { return parseInt(hex.slice(i, i + 2), 16); });
    return 'rgba(' + v.join(',') + ',' + a + ')';
  }

  // ================================================================ plotly helpers
  var registry = [];
  var CONFIG = { displaylogo: false, responsive: true, modeBarButtonsToRemove: ['autoScale2d', 'lasso2d', 'select2d'] };

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
      legend: Object.assign({ orientation: 'h', y: -0.2, font: { color: c.secondary } }, o.legendOpts || {}),
      margin: o.margin || { l: 56, r: 16, t: 16, b: 48 },
      hoverlabel: { font: { family: FONT } },
      hovermode: o.hovermode || 'closest',
      annotations: o.annotations, shapes: o.shapes, barmode: o.barmode, dragmode: o.dragmode, bargap: o.bargap
    };
  }
  function plot(div, build) {
    if (div._note) { clear(div); div._note = false; }
    div._render = function () {
      var s = build(resolved());
      var cfg = Object.assign({}, CONFIG, s.config || {});
      s.layout.autosize = true;
      window.Plotly.react(div, s.data, s.layout, cfg);
      if (s.onClick && div.on) {
        if (div.removeAllListeners) div.removeAllListeners('plotly_click');
        div.on('plotly_click', s.onClick);
      }
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
    if (window.Plotly && div.data) window.Plotly.purge(div);
    clear(div); div._render = null; div._note = true; div._dirty = false;
    div.appendChild(emptyNote(text));
  }
  function refreshVisible(root) {
    $$('.plot', root).forEach(function (d) {
      if (!d._render) return;
      if (d._dirty) d._render(); else if (window.Plotly && d.data) window.Plotly.Plots.resize(d);
    });
  }
  function plotBox(cls) { return h('div', { class: 'plot ' + (cls || '') }); }
  function emptyNote(text) { return h('div', { class: 'empty', text: text }); }

  // ================================================================ widgets
  function chip(level, text) {
    return h('span', { class: 'chip ' + level },
      h('span', { class: 'ico', 'aria-hidden': 'true', text: ICON[level] || '' }), text === undefined ? LEVEL_LABEL[level] : text);
  }
  function confChip(v) { return h('span', { class: 'chip conf ' + (v || 'NA') }, v === 'NA' || !v ? 'n/a' : v); }
  function typeChip(v) {
    if (v !== 'DEL' && v !== 'DUP') return fmt(v);
    return h('span', { class: 'chip type ' + v }, h('span', { class: 'ico', 'aria-hidden': 'true', text: v === 'DEL' ? '\u25BC' : '\u25B2' }), v);
  }
  function sevChip(v) { return v ? chip(v, v) : ''; }
  function card(title, hint, body, cls) {
    return h('section', { class: 'card' + (cls ? ' ' + cls : '') },
      h('header', null, h('h3', { text: title }), hint ? h('span', { class: 'muted small hint', text: hint }) : null), body);
  }
  function detailCard(title, hint, body, open) {
    var details = h('details', { class: 'card collapsible-section' });
    var summary = h('summary', null, h('span', { text: title }), hint ? h('span', { class: 'collapse-hint', text: hint }) : null);
    details.appendChild(summary);
    details.appendChild(h('div', { class: 'collapsible-body' }, body));
    if (open) details.open = true;
    return details;
  }
  function tile(value, label, level, detail) {
    return h('div', { class: 'tile' + (level ? ' ' + level : '') }, h('div', { class: 'v', text: String(value) }),
      h('div', { class: 'l', text: label }), detail ? h('span', { class: 'detail', text: detail }) : null);
  }
  function alertBox(level, title, text, extra) {
    return h('div', { class: 'alert ' + level, role: level === 'fail' || level === 'error' ? 'alert' : 'note' },
      h('span', { class: 'ico', 'aria-hidden': 'true', text: ICON[level] }),
      h('div', { class: 'body' },
        title ? h('span', { class: 'ttl', text: title }) : null, title && text ? ' ' : null, text || null, extra || null));
  }
  function kv(pairs) {
    var dl = h('dl', { class: 'kv' });
    pairs.forEach(function (p) {
      if (p[1] === null || p[1] === undefined || p[1] === '') return;
      dl.appendChild(h('dt', { text: p[0] })); dl.appendChild(h('dd', null, p[1]));
    });
    return dl;
  }
  function runKv() {
    return kv([
      ['Report version', META.version],
      ['Genome build', META.genome_build || META.genome_version],
      ['Generated', META.generated]
    ]);
  }
  function field(label, el) { return h('label', { class: 'field' }, label, el); }
  function selectOf(options, value, onchange, aria) {
    var sel = h('select', { 'aria-label': aria || 'Select' }, options.map(function (o) {
      return typeof o === 'string' ? h('option', { value: o, text: o }) : h('option', { value: o.value, text: o.text });
    }));
    if (value !== undefined) sel.value = value;
    if (onchange) sel.addEventListener('change', function () { onchange(sel.value); });
    return sel;
  }

  // ================================================================ tables
  function colType(T, c) {
    var arr = tcol(T, c);
    for (var i = 0; i < arr.length; i++) { if (arr[i] !== null && arr[i] !== undefined) return typeof arr[i]; }
    return 'string';
  }
  // T = {cols, data: {col: [..]}, n}. opts: filename, pageSize, render {col: fn(v,i)->node},
  // cellClass {col: fn(v,i)->class}, labels {col: header}, sort [col,'asc'|'desc'],
  // onRow fn(i), selected fn()->row index, wrap [cols], note
  function tableView(host, T, opts) {
    opts = opts || {};
    var n = T.n, pageSize = opts.pageSize || 25, page = 0, filter = '', pred = null;
    var cols = (T.cols || []).slice();
    var sortCol = opts.sort ? opts.sort[0] : null, sortDir = opts.sort ? opts.sort[1] : 'asc';
    var types = {}; cols.forEach(function (c) { types[c] = colType(T, c); });
    var idx = range(n);
    var search = h('input', { type: 'search', placeholder: 'Filter rows\u2026', 'aria-label': 'Filter rows' });
    var count = h('span', { class: 'muted small' });
    var pager = h('span', { class: 'pager' });
    var tbody = h('tbody'), thead = h('thead');
    var table = h('table', { class: 'tbl' }, thead, tbody);
    var wrap = h('div', { class: 'tbl-wrap' }, table);
    var dl = h('button', { class: 'btn', type: 'button', text: 'Download CSV', onclick: function () {
      download((opts.filename || 'table') + '.csv', toCSV(cols, function (c, k) { return tcol(T, c)[idx[k]]; }, idx.length));
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
          for (var j = 0; j < cols.length && !hit; j++) {
            var v = tcol(T, cols[j])[i];
            if (v !== null && v !== undefined && String(v).toLowerCase().indexOf(f) >= 0) hit = true;
          }
          if (!hit) continue;
        }
        idx.push(i);
      }
      if (sortCol) {
        var arr = tcol(T, sortCol), dir = sortDir === 'asc' ? 1 : -1, num = types[sortCol] === 'number';
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
      cols.forEach(function (c) {
        tr.appendChild(h('th', {
          scope: 'col', tabindex: '0', text: (opts.labels && opts.labels[c]) || pretty(c),
          'aria-sort': sortCol === c ? (sortDir === 'asc' ? 'ascending' : 'descending') : 'none',
          onclick: function () { sortDir = sortCol === c && sortDir === 'asc' ? 'desc' : 'asc'; sortCol = c; recompute(); },
          onkeydown: function (e) { if (e.key === 'Enter') this.click(); }
        }));
      });
      clear(thead).appendChild(tr);
      clear(tbody);
      var sel = opts.selected ? opts.selected() : -1;
      var start = page * pageSize, end = Math.min(idx.length, start + pageSize);
      var wrapCols = opts.wrap || [];
      for (var k = start; k < end; k++) {
        var i = idx[k];
        var row = h('tr', { class: (opts.onRow ? 'clickable ' : '') + (i === sel ? 'sel' : '') });
        if (opts.onRow) (function (ii) {
          row.tabIndex = 0;
          row.addEventListener('click', function (e) { if (e.target.closest && e.target.closest('button,a')) return; opts.onRow(ii); });
          row.addEventListener('keydown', function (e) { if (e.key === 'Enter') opts.onRow(ii); });
        })(i);
        cols.forEach(function (c) {
          var v = tcol(T, c)[i], r = opts.render && opts.render[c];
          var cc = opts.cellClass && opts.cellClass[c] ? opts.cellClass[c](v, i) : '';
          var cls = ((types[c] === 'number' ? 'num ' : '') + (wrapCols.indexOf(c) >= 0 ? 'wrap ' : '') + (cc || '')).trim();
          row.appendChild(h('td', { class: cls || false }, r ? r(v, i) : fmt(v)));
        });
        tbody.appendChild(row);
      }
      count.textContent = fmtInt(idx.length) + (idx.length === n ? '' : ' of ' + fmtInt(n)) + ' rows' +
        (T.truncated ? ' (truncated from ' + fmtInt(T.n_total) + ')' : '');
      var pages = Math.max(1, Math.ceil(idx.length / pageSize));
      clear(pager);
      pager.appendChild(h('button', { class: 'btn', type: 'button', text: '\u2039', 'aria-label': 'Previous page', disabled: page === 0, onclick: function () { page--; draw(); } }));
      pager.appendChild(h('span', { class: 'small', text: (page + 1) + ' / ' + pages }));
      pager.appendChild(h('button', { class: 'btn', type: 'button', text: '\u203a', 'aria-label': 'Next page', disabled: page >= pages - 1, onclick: function () { page++; draw(); } }));
    }
    var timer = null;
    search.addEventListener('input', function () {
      clearTimeout(timer);
      timer = setTimeout(function () { filter = search.value; page = 0; recompute(); }, 120);
    });
    recompute();
    return { setPredicate: function (fn) { pred = fn; page = 0; recompute(); }, refresh: recompute, redraw: draw, count: function () { return idx.length; }, visible: function () { return idx.slice(); }, select: function (i) { draw(); } };
  }
  function simpleTable(T, opts) { var host = h('div'); tableView(host, T, opts); return host; }
  function subTable(T, cols) {
    var out = { cols: cols, data: {}, n: T.n };
    cols.forEach(function (c) { out.data[c] = tcol(T, c); });
    return out;
  }

  // ================================================================ global state & filters
  var state = { sample: META.focus || '', call: null, alertFilter: 'all', listeners: [] };
  var headerPicker = null, globalSeverityControl = null;
  function onState(cb) { state.listeners.push(cb); }
  function severityGroup(level) {
    if (level === 'error' || level === 'fail') return 'fail';
    if (level === 'warn' || level === 'warning') return 'warning';
    if (level === 'info') return 'note';
    if (level === 'ok') return 'ok';
    return 'note';
  }
  function severityMatches(level, filter) {
    if (!filter || filter === 'all') return true;
    if (filter === 'hide') return false;
    return severityGroup(level) === filter;
  }
  function applyGlobalAlertFilter() {
    var filter = state.alertFilter || 'all';
    document.body.setAttribute('data-alert-filter', filter);
    $$('.alert').forEach(function (el) {
      var classes = (el.className || '').split(/\s+/);
      var level = classes.filter(function (c) { return /^(fail|error|warn|warning|info|ok)$/.test(c); })[0] || 'info';
      el.hidden = !severityMatches(level, filter);
    });
    if (globalSeverityControl) globalSeverityControl.set(filter);
  }
  function notify() {
    state.listeners.forEach(function (cb) { try { cb(); } catch (e) { if (window.console) console.error(e); } });
    applyGlobalAlertFilter();
  }
  // A single cohort-wide sample filter lives in the report header.
  function setSample(s, quiet) {
    var next = s || '';
    if (next !== state.sample) state.call = null;
    state.sample = next;
    if (headerPicker) headerPicker.value = state.sample;
    if (!quiet) notify();
  }
  function setAlertFilter(filter) {
    state.alertFilter = filter || 'all';
    notify();
  }
  function makeSeveritySegments(current, options, onChange, ariaLabel) {
    var buttons = {}, seg = h('div', { class: 'seg', role: 'group', 'aria-label': ariaLabel || 'Filter by severity' });
    options.forEach(function (opt) {
      var b = h('button', { type: 'button', 'aria-pressed': opt.value === current ? 'true' : 'false', text: opt.label,
        onclick: function () { current = opt.value; sync(); if (onChange) onChange(current); } });
      buttons[opt.value] = b; seg.appendChild(b);
    });
    function sync() { Object.keys(buttons).forEach(function (k) { buttons[k].setAttribute('aria-pressed', k === current ? 'true' : 'false'); }); }
    return { el: seg, set: function (value) { current = value; sync(); } };
  }

  // ================================================================ sub-tabs
  // Builds an inner segmented control + panels.  items: [{id, label, build(root)}]
  function subTabs(items, defaultId) {
    var current = defaultId || items[0].id;
    var bar = h('div', { class: 'subtabs', role: 'tablist' });
    var panels = {};
    var container = h('div');
    items.forEach(function (it) {
      var btn = h('button', { class: 'subtab', role: 'tab', type: 'button', 'aria-selected': it.id === current ? 'true' : 'false',
        text: it.label, onclick: function () { switchTo(it.id); } });
      bar.appendChild(btn);
      var p = h('div', { class: 'subpanel', role: 'tabpanel', hidden: it.id === current ? false : true });
      panels[it.id] = { btn: btn, panel: p, built: false, build: it.build };
      container.appendChild(p);
    });
    function switchTo(id) {
      current = id;
      Object.keys(panels).forEach(function (k) {
        var on = k === id; panels[k].panel.hidden = !on;
        panels[k].btn.setAttribute('aria-selected', on ? 'true' : 'false');
      });
      if (!panels[id].built) { panels[id].built = true; try { panels[id].build(panels[id].panel); } catch (e) { panels[id].panel.appendChild(h('div', { class: 'empty', text: 'Error: ' + e.message })); if (window.console) console.error(e); } }
      else refreshVisible(panels[id].panel);
      applyGlobalAlertFilter();
      setTimeout(function () { refreshVisible(panels[id].panel); applyGlobalAlertFilter(); }, 60);
    }
    container.insertBefore(bar, container.firstChild);
    // build the default panel immediately
    if (!panels[current].built) { panels[current].built = true; panels[current].build(panels[current].panel); }
    return { el: container, switchTo: switchTo, refresh: function () { if (panels[current].built) refreshVisible(panels[current].panel); } };
  }

  // ================================================================ call helpers
  var CALL_ROW = {};
  tcol(CALLS, 'id').forEach(function (id, i) { CALL_ROW[id] = i; });
  if (!hasCol(CALLS, 'id') && hasCol(CALLS, 'ID')) tcol(CALLS, 'ID').forEach(function (id, i) { CALL_ROW[id] = i; });

  function callsOf(s) {
    var out = [], sm = tcol(CALLS, 'sample');
    if (!sm.length) sm = tcol(CALLS, 'Sample');
    for (var i = 0; i < CALLS.n; i++) if (!s || sm[i] === s) out.push(i);
    return out;
  }
  var CONF_RANK = { HIGH: 3, MEDIUM: 2, LOW: 1 };
  function callOrder(rows) {
    var q = tcol(CALLS, 'q_some'), bf = tcol(CALLS, 'bf'), cf = tcol(CALLS, 'confidence');
    if (!cf.length) cf = tcol(CALLS, 'Confidence');
    return rows.slice().sort(function (a, b) {
      var d = (CONF_RANK[cf[b]] || 0) - (CONF_RANK[cf[a]] || 0);
      if (d) return d;
      var qa = isNum(q[a]) ? q[a] : isNum(bf[a]) ? bf[a] : -1;
      var qb = isNum(q[b]) ? q[b] : isNum(bf[b]) ? bf[b] : -1;
      return qb - qa;
    });
  }
  function hasPlot(i) {
    var hp = tcol(CALLS, 'has_plot'); if (hp.length) return !!hp[i];
    var v = tcol(CALLS, 'view'); return !!v[i];
  }
  function callId(i) { return tcol(CALLS, 'id')[i] !== undefined ? tcol(CALLS, 'id')[i] : tcol(CALLS, 'ID')[i]; }
  function callSample(i) { return tcol(CALLS, 'sample')[i] || tcol(CALLS, 'Sample')[i]; }
  function callType(i) { return tcol(CALLS, 'type')[i] || tcol(CALLS, 'Type')[i]; }
  function callConf(i) { return tcol(CALLS, 'confidence')[i] || tcol(CALLS, 'Confidence')[i]; }
  function callGene(i) { return tcol(CALLS, 'gene')[i] || tcol(CALLS, 'Gene')[i] || ''; }
  function callInterval(i) { return tcol(CALLS, 'interval')[i] || tcol(CALLS, 'Genomic')[i] || ''; }

  // ================================================================ TABS registry
  // Keep the tab registry separate from call-selection helpers: the shell
  // uses it to create and activate the report panels during init().
  var TABS = [];
  function tab(id, label, available, build) {
    var t = { id: id, label: label, available: available, build: build, built: false };
    TABS.push(t); return t;
  }

  // ================================================================ OVERVIEW
  function overviewBars() {
    var host = h('div'), box = plotBox('');
    var order = selectOf([{ value: 'input', text: 'Input order' }, { value: 'calls', text: 'Most calls first' }, { value: 'flags', text: 'Most flags first' }], 'calls', null, 'Sample order');
    host.appendChild(h('div', { class: 'toolbar' }, field('Order', order)));
    host.appendChild(box);
    if (!SAMPLES.length) { host.appendChild(emptyNote('No samples.')); return host; }
    function draw() {
      var nC = tcol(S, 'n_cnvs'), nC2 = hasCol(S, 'CNV calls') ? tcol(S, 'CNV calls') : nC;
      var nH = tcol(S, 'n_high'), nH2 = hasCol(S, 'HIGH') ? tcol(S, 'HIGH') : nH;
      var nM = tcol(S, 'n_medium'), nM2 = hasCol(S, 'MEDIUM') ? tcol(S, 'MEDIUM') : nM;
      var nL = tcol(S, 'n_low'), nL2 = hasCol(S, 'LOW') ? tcol(S, 'LOW') : nL;
      var ids = range(S.n);
      if (order.value === 'calls') ids.sort(function (a, b) { return (nC2[b] || 0) - (nC2[a] || 0); });
      plot(box, function (t) {
        var c = SEM[t], x = ids.map(function (i) { return SAMPLES[i]; });
        var mk = function (name, arr, col) { return { type: 'bar', name: name, x: x, y: ids.map(function (i) { return arr[i] || 0; }), marker: { color: col }, hovertemplate: '%{x}<br>' + name + ': %{y}<extra></extra>' }; };
        var data = [mk('HIGH', nH2.length ? nH2 : nH, c.HIGH || c.high), mk('MEDIUM', nM2.length ? nM2 : nM, c.MEDIUM || c.medium), mk('LOW', nL2.length ? nL2 : nL, c.LOW || c.low)];
        var L = baseLayout(t, { ylab: 'CNV calls', barmode: 'stack', legend: true, xaxis: { type: 'category', tickangle: -45 }, margin: { l: 56, r: 16, t: 16, b: 90 } });
        L.yaxis.rangemode = 'tozero';
        return { data: data, layout: L, onClick: function (ev) { var p = ev.points && ev.points[0]; if (p) setSample(p.x); } };
      });
    }
    order.addEventListener('change', draw); draw();
    return host;
  }
  function healthStatusTiles() {
    var counts = { ok: 0, warn: 0, fail: 0 };
    var lvCol = tcol(S, 'level'), stCol = tcol(S, 'status');
    if (lvCol.length) lvCol.forEach(function (l) {
      var key = String(l || '').toLowerCase();
      if (key === 'warning' || key === 'excluded') key = 'warn';
      if (key === 'error') key = 'fail';
      if (counts[key] !== undefined) counts[key]++;
    });
    else if (stCol.length) stCol.forEach(function (l) {
      var key = ({ OK: 'ok', Warning: 'warn', Fail: 'fail', Excluded: 'warn', failed: 'fail', excluded: 'warn' })[l] || 'ok';
      counts[key]++;
    });
    var grid = h('div', { class: 'health-tile-grid' });
    [
      ['ok', 'No flags', 'Samples without QC flags'],
      ['warn', 'Warnings', 'Samples with QC warnings'],
      ['fail', 'Fail', 'Samples failing QC']
    ].forEach(function (p) {
      var pct = S.n ? Math.round(100 * counts[p[0]] / S.n) : 0;
      var el = tile(counts[p[0]], p[1], p[0], pct + '% of samples');
      el.classList.add('health-tile');
      el.title = p[2];
      grid.appendChild(el);
    });
    return grid;
  }

  function missingBedChromosomeClasses() {
    var missing = [];
    if (SCOPE && SCOPE.n) {
      range(SCOPE.n).forEach(function (i) {
        var cls = String(tcol(SCOPE, 'class')[i] || '').toLowerCase();
        var present = !!tcol(SCOPE, 'present')[i];
        if (!present && cls !== 'other' && (cls === 'autosome' || cls === 'x' || cls === 'y')) {
          var label = cls === 'autosome' ? 'Autosomes (chr1-22)' : 'chr' + cls.toUpperCase();
          missing.push(label);
        }
      });
    } else if (CHR.status && CHR.status.n) {
      range(CHR.status.n).forEach(function (i) {
        var cls = String(tcol(CHR.status, 'class')[i] || '').toLowerCase();
        var inInput = !!tcol(CHR.status, 'in_input')[i];
        if (!inInput && cls !== 'other' && (cls === 'autosome' || cls === 'x' || cls === 'y')) {
          missing.push(cls === 'autosome' ? 'Autosomes (chr1-22)' : 'chr' + cls.toUpperCase());
        }
      });
    }
    return uniq(missing);
  }

  var overviewSubtabs = null;
  function buildOverviewSamples(root) {
    var del = 0, dup = 0, high = 0;
    var typeCol = tcol(CALLS, 'type'); if (!typeCol.length) typeCol = tcol(CALLS, 'Type');
    typeCol.forEach(function (t) { if (t === 'DEL') del++; else if (t === 'DUP') dup++; });
    var confCol = tcol(CALLS, 'confidence'); if (!confCol.length) confCol = tcol(CALLS, 'Confidence');
    confCol.forEach(function (c) { if (c === 'HIGH') high++; });
    var called = Number(SUM.samples_called !== undefined ? SUM.samples_called : (META.n_analysed || 0));
    var requested = Number(SUM.samples_requested !== undefined ? SUM.samples_requested : (META.n_samples || S.n));
    var notAnalysed = Math.max(0, requested - called);
    // Overview > Samples contains the cohort snapshot and the sample QC table.
    // Cohort-wide Sample QC status tiles remain in QC & warnings.
    var summary = h('div', { class: 'overview-summary overview-summary-single' });
    var stats = h('div', { class: 'overview-summary-stats' },
      h('div', { class: 'tiles' },
        tile(called + ' / ' + requested, 'samples analysed / requested', called < requested ? 'fail' : 'ok'),
        tile(notAnalysed, 'samples not analysed', notAnalysed > 0 ? 'fail' : 'ok'),
        tile(fmtInt(CALLS.n), 'CNV calls (' + del + ' DEL, ' + dup + ' DUP)'),
        tile(fmtInt(high), 'HIGH-confidence calls'),
        tile(SUM.samples_with_calls || 0, 'samples with at least one call'),
        tile(fmtInt(META.n_targets || SUM.targets || 0), 'targets analysed')));
    // Resolve the names from sample-level analysis/status fields.
    var sampleData = S.data || {}, sampleNames = sampleData.sample || sampleData.Sample || SAMPLES;
    var analysed = sampleData.Analysed || sampleData.analysed || [];
    var statuses = sampleData.status || [];
    var notAnalysedNames = sampleNames.filter(function (name, i) {
      return (analysed.length && /^(no|false|0)$/i.test(String(analysed[i]))) ||
        (!analysed.length && /^(excluded|failed|not analysed|not analyzed)$/i.test(String(statuses[i] || '')));
    });
    summary.appendChild(stats);
    root.appendChild(card('Cohort snapshot', 'Cohort-level summary and analysis completeness', summary));
    buildSampleTable(root);
  }
  function buildOverviewCalls(root) {
    root.appendChild(card('Calls per sample', 'Select a bar to filter the report to that sample', overviewBars()));
    if (CALLS.n) {
      var overviewMap = plotBox('');
      root.appendChild(card('Genome-wide map of calls', 'DEL and DUP calls per chromosome; filtered by the selected sample', overviewMap));
      genomeCallHistogram(overviewMap);
      onState(function () { if (overviewMap._render) overviewMap._render(); });
    } else root.appendChild(emptyNote('No CNV calls are available for the genome-wide distribution.'));
  }
  tab('overview', 'Overview', function () { return true; }, function (root) {
    var st = subTabs([
      { id: 'samples', label: 'Samples', build: buildOverviewSamples },
      { id: 'chromosomes', label: 'Chromosomes', build: buildChromosomeOverview },
      { id: 'calls', label: 'Call landscape', build: buildOverviewCalls }
    ], 'samples');
    overviewSubtabs = st;
    st.el.firstChild.setAttribute('aria-label', 'Overview sections');
    root.appendChild(st.el);
    onState(function () { st.refresh(); });
  });

  // ================================================================ QC & WARNINGS
  // QC_METRICS for the metric bars
  var QC_METRICS = [
    { col: 'median_depth', col2: 'Median depth', label: 'Median target depth', thr: function () { return TH.qc_min_cov || TH.min_cov; }, dir: 'below', hint: 'below the line = low coverage' },
    { col: 'total_reads', col2: 'Total reads', label: 'Total reads on the analysed targets', thr: function () { return TH.qc_min_total_reads || TH.min_total_reads; }, dir: 'below', hint: 'below the line = too few reads' },
    { col: 'max_corr', col2: 'Max correlation', label: 'Best Spearman correlation with another sample', thr: function () { return TH.qc_min_corr || TH.min_corr; }, dir: 'below', hint: 'below the line = no comparable reference' },
    { col: 'noise_z', col2: '', label: 'Sample noise (robust z of the exon-count dispersion)', thr: function () { return TH.qc_zscore; }, dir: 'above', hint: 'above the line = noisy sample' },
    { col: 'n_refs', col2: 'References', label: 'Reference samples used', thr: function () { return IS_CANOPE ? 5 : null; }, dir: 'below', hint: 'below the line = weak reference set' },
    { col: 'n_cnvs', col2: 'CNV calls', label: 'CNV calls', thr: function () { return null; }, dir: null, hint: 'outliers suggest a noisy sample / batch effect' }
  ];
  function qcTypeChip(v) {
    var lv = v === 'Missing Chromosome' ? 'fail' : v === 'Chromosome not analysed' ? 'warn' : v === 'Whole sample' ? 'warn' : 'info';
    return h('span', { class: 'chip ' + lv }, h('span', { class: 'ico', 'aria-hidden': 'true', text: ICON[lv] }), v);
  }
  function metricBars(m, orderSel) {
    var box = plotBox('short');
    var vals = tcol(S, m.col);
    if (!vals.length && m.col2) vals = tcol(S, m.col2);
    var ids = range(S.n).filter(function (i) { return isNum(vals[i]); });
    if (!ids.length) return null;
    function draw() {
      var o = orderSel.value, order = ids.slice();
      if (o === 'asc') order.sort(function (a, b) { return vals[a] - vals[b]; });
      else if (o === 'desc') order.sort(function (a, b) { return vals[b] - vals[a]; });
      plot(box, function (t) {
        var c = SEM[t], thr = m.thr();
        var bad = function (i) { return thr !== null && thr !== undefined && (m.dir === 'below' ? vals[i] < thr : vals[i] > thr); };
        var shapes = (thr !== null && thr !== undefined) ? [{ type: 'line', xref: 'paper', yref: 'y', x0: 0, x1: 1, y0: thr, y1: thr, line: { color: c.fail, dash: 'dash', width: 1.5 } }] : [];
        var d = [{ type: 'bar', x: order.map(function (i) { return SAMPLES[i]; }), y: order.map(function (i) { return vals[i]; }),
          marker: { color: order.map(function (i) { return bad(i) ? c.warn : CAT[t][0]; }) },
          hovertemplate: '%{x}<br>' + m.label + ': %{y}<extra></extra>' }];
        var lo = Math.min(0, Math.min.apply(null, ids.map(function (i) { return vals[i]; }))), hi = Math.max(Math.max.apply(null, ids.map(function (i) { return vals[i]; })), thr !== null && thr !== undefined ? thr : -Infinity);
        var L = baseLayout(t, { ylab: '', yaxis: { type: 'linear', range: [lo < 0 ? lo * 1.1 : 0, hi > 0 ? hi * 1.08 : 1] }, xaxis: { type: 'category', tickangle: -45, showticklabels: order.length <= 60 }, shapes: shapes, margin: { l: 60, r: 16, t: 10, b: order.length <= 60 ? 80 : 20 } });
        return { data: d, layout: L, onClick: function (ev) { var p = ev.points && ev.points[0]; if (p) setSample(p.x); } };
      });
    }
    orderSel.addEventListener('change', draw);
    draw();
    var nbad = m.thr() === null || m.thr() === undefined ? 0 : ids.filter(function (i) { return m.dir === 'below' ? vals[i] < m.thr() : vals[i] > m.thr(); }).length;
    return card(m.label, (m.thr() !== null && m.thr() !== undefined ? 'threshold ' + fmt(m.thr()) + ' \u2013 ' : '') + m.hint + (nbad ? ' \u2013 ' + nbad + ' sample(s) flagged' : ''), box, nbad ? 'flagged' : '');
  }

  // ---- Cohort metrics ----
  function buildCohortMetricsSubtab(root, metrics) {
    var stack = h('div', { class: 'stack' });
    var orderSel = selectOf([{ value: 'input', text: 'Input order' }, { value: 'asc', text: 'Lowest first' }, { value: 'desc', text: 'Highest first' }], 'input', null, 'Order');
    stack.appendChild(h('div', { class: 'toolbar' }, field('Sort bars', orderSel), h('span', { class: 'note', text: 'Dashed lines mark QC thresholds. Select a bar to inspect that sample across the report.' })));
    var grid = h('div', { class: 'grid g2' });
    metrics.forEach(function (m) { var c = metricBars(m, orderSel); if (c) grid.appendChild(c); });
    if (grid.children.length) stack.appendChild(grid);
    else stack.appendChild(emptyNote('No metrics in this group were provided.'));
    root.appendChild(stack);
  }
  function buildPcaSubtab(root) {
    if (!PCA) { root.appendChild(emptyNote('PCA needs at least 3 samples with variable coverage.')); return; }
    var pbox = plotBox('tall');
    root.appendChild(card('Cohort PCA of coverage profiles', (PCA.n_samples || PCA.names && PCA.names.length) + ' samples · ' + fmtInt(PCA.n_targets || PCA.n_rows || 0) + ' targets (log2 counts)', pbox));
    pcaPlot(pbox);
  }
  function buildDepthHistogramSubtab(root) {
    if (!DEPTH_HIST || !DEPTH_HIST.counts || !DEPTH_HIST.counts.length) {
      root.appendChild(emptyNote('No target-depth distribution was provided.')); return;
    }
    var dbox = plotBox('mid');
    root.appendChild(card('Target depth distribution', 'Median depth of each target; dashed line = QC threshold', dbox));
    depthHist(dbox);
  }

  function buildDepthCorrelationSubtab(root) {
    var stack = h('div', { class: 'stack' });
    var depth = tcol(S, 'median_depth'); if (!depth.length) depth = tcol(S, 'Median depth');
    var corr = tcol(S, 'max_corr'); if (!corr.length) corr = tcol(S, 'Max correlation');
    var okIds = range(S.n).filter(function (i) { return isNum(depth[i]) && isNum(corr[i]); });
    if (okIds.length >= 2) {
      var sc = plotBox('tall');
      stack.appendChild(card('Depth vs. best correlation', 'dashed lines = QC thresholds; click a point to select the sample', sc));
      plot(sc, function (t) {
        var c = SEM[t], data = [];
        var lvCol = tcol(S, 'level'), stCol = tcol(S, 'status');
        ['ok', 'warn', 'fail'].forEach(function (lv) {
          var ids = okIds.filter(function (i) {
            if (lvCol.length) return lvCol[i] === lv;
            if (stCol.length) { var m = { OK: 'ok', Warning: 'warn', Fail: 'fail', Excluded: 'warn' }[stCol[i]] || 'ok'; return m === lv; }
            return lv === 'ok';
          });
          if (!ids.length) return;
          data.push({ type: 'scatter', mode: 'markers', name: LEVEL_LABEL[lv], x: ids.map(function (i) { return depth[i]; }), y: ids.map(function (i) { return corr[i]; }),
            text: ids.map(function (i) { return SAMPLES[i]; }), customdata: ids.map(function (i) { return SAMPLES[i]; }),
            marker: { color: c[lv], size: 9, symbol: lv === 'ok' ? 'circle' : lv === 'warn' ? 'diamond' : 'x', line: { width: 1, color: CHROME[t].surface } },
            hovertemplate: '<b>%{text}</b><br>depth %{x}<br>corr %{y}<extra>' + LEVEL_LABEL[lv] + '</extra>' });
        });
        var shapes = [];
        if (TH.qc_min_cov || TH.min_cov) shapes.push({ type: 'line', xref: 'x', yref: 'paper', x0: TH.qc_min_cov || TH.min_cov, x1: TH.qc_min_cov || TH.min_cov, y0: 0, y1: 1, line: { color: c.fail, dash: 'dash' } });
        if (TH.qc_min_corr || TH.min_corr) shapes.push({ type: 'line', xref: 'paper', yref: 'y', x0: 0, x1: 1, y0: TH.qc_min_corr || TH.min_corr, y1: TH.qc_min_corr || TH.min_corr, line: { color: c.fail, dash: 'dash' } });
        return { data: data, layout: baseLayout(t, { xlab: 'median depth', ylab: 'best correlation', legend: true, shapes: shapes }), onClick: function (ev) { var p = ev.points && ev.points[0]; if (p && p.customdata) setSample(p.customdata); } };
      });
    } else stack.appendChild(emptyNote('At least two samples with both depth and correlation values are needed.'));
    root.appendChild(stack);
  }

  // ---- sub-tab: Flags ----
  function buildSampleQcFlags(root) {
    var stack = h('div', { class: 'stack' });
    // ECHO QC flags table
    if (QC_FLAGS && QC_FLAGS.n) {
      var host = h('div');
      var typ = '';
      var tv = null;
      var bar = h('div', { class: 'toolbar' });
      var types = uniq(tcol(QC_FLAGS, 'Type')).sort();
      var typeSel = h('select', { 'aria-label': 'Flag type' }, h('option', { value: '', text: 'All types' }), types.map(function (x) { return h('option', { value: x, text: x }); }));
      typeSel.addEventListener('change', function () { typ = typeSel.value; apply(); });
      bar.appendChild(h('span', { class: 'note', text: 'This table is filtered by sample and flag type; the global severity selector does not affect it.' }));
      bar.appendChild(field('Flag type', typeSel));
      host.appendChild(bar);
      var tblHost = h('div'); host.appendChild(tblHost);
      tv = tableView(tblHost, QC_FLAGS, { filename: 'qc_flags', pageSize: 20, wrap: ['Details'], sort: ['Severity', 'asc'],
        render: { Severity: function (v) { return sevChip(v); }, Details: function (v) { return v; } } });
      function apply() {
        var bySample = function (i) { var s = state.sample; if (!s) return true; var v = tcol(QC_FLAGS, 'Sample')[i]; return v === s || v === 'All'; };
        tv.setPredicate(function (i) {
          return bySample(i) && (!typ || tcol(QC_FLAGS, 'Type')[i] === typ);
        });
      }
      apply();
      // React only to sample changes. Severity changes must not alter this table.
      var lastSample = state.sample;
      onState(function () {
        if (lastSample !== state.sample) { lastSample = state.sample; apply(); }
      });
      stack.appendChild(card('Sample QC flags', 'Whole sample \u00b7 Whole exon \u00b7 Exon variability \u00b7 Cohort size \u00b7 Missing Chromosome', host));
    }
    // CANOPE QC metrics block
    if (QC && QC.rows && QC.rows.n) {
      var host2 = h('div'), typeBox = plotBox('short');
      var tcounts = tcol(QC.types, 'Freq'), tnames = tcol(QC.types, 'Type');
      stack.appendChild(h('div', { class: 'grid g-main' },
        card('Exon-level QC flags by type', QC.n_total + ' flags in total', tnames.length ? typeBox : emptyNote('No QC flags.')),
        card('Genes with the most QC flags', 'top 60', QC.genes && QC.genes.n ? simpleTable(QC.genes, { filename: 'qc_genes', pageSize: 10 }) : emptyNote('No exon-level flags.'))));
      if (tnames.length) plot(typeBox, function (t) {
        return { data: [{ type: 'bar', orientation: 'h', y: tnames, x: tcounts, marker: { color: CAT[t][3] }, hovertemplate: '%{y}: %{x}<extra></extra>' }],
          layout: baseLayout(t, { xlab: 'flags', yaxis: { autorange: 'reversed' }, margin: { l: 170, r: 16, t: 10, b: 40 } }) };
      });
      var tv2 = tableView(host2, QC.rows, { filename: 'qc_metrics_flags', pageSize: 20, render: { Type: qcTypeChip },
        note: QC.capped ? 'Table capped at ' + fmtInt(QC.rows.n) + ' rows; full list in the QC TSV.' : null });
      var sm = tcol(QC.rows, 'Sample');
      var apply2 = function () { tv2.setPredicate(state.sample ? function (i) { return sm[i] === state.sample || sm[i] === 'All'; } : null); };
      onState(apply2); apply2();
      stack.appendChild(card('Exon-level QC details', 'Filtered by the selected sample', host2));
    }
    if (!stack.children.length) stack.appendChild(card('Sample QC flags', null, emptyNote('No QC metrics file was provided.')));
    root.appendChild(stack);
  }

  // ---- Chromosome scope: one canonical summary, preferring target counts when available ----
  function chromosomeScopeCard() {
    if (SCOPE && SCOPE.n) {
      var scopeTiles = range(SCOPE.n).map(function (i) {
        var present = !!tcol(SCOPE, 'present')[i], analysed = !!tcol(SCOPE, 'analysed')[i];
        var status = !present ? 'Not in BED' : analysed ? 'Analysed' : 'In BED, not analysed';
        var chrClass = String(tcol(SCOPE, 'class')[i] || '').toLowerCase();
        var cls = !present && chrClass !== 'other' ? 'fail' : !present ? '' : analysed ? 'ok' : 'fail';
        return tile(fmtInt(tcol(SCOPE, 'targets')[i]), tcol(SCOPE, 'label')[i], cls, status);
      });
      return card('Run scope by chromosome class', META.mode_label || META.modechrom_label || ('Mode ' + (META.modechrom || 'n/a')),
        h('div', null, h('div', { class: 'scope' }, scopeTiles),
          h('p', { class: 'note', text: 'A run analyses one chromosome mode. Classes present in the BED but not analysed require a separate run; absence of calls is not a negative result for those regions.' })));
    }
    if (CHR.status && CHR.status.n) {
      return card('Run scope by chromosome class', META.mode_label || META.modechrom_label, chromStrip());
    }
    return null;
  }

  tab('qc', 'QC & warnings', function () { return true; }, function (root) {
    // Sample QC status belongs with the severity/warning view, not Overview.
    root.appendChild(card('Sample QC status', 'Cohort-wide distribution of sample QC levels', healthStatusTiles()));
    // Sample-level QC flags remain independent from the global severity filter.
    buildSampleQcFlags(root);
    var st = subTabs([
      { id: 'coverage', label: 'Coverage', build: function (r) { buildCohortMetricsSubtab(r, QC_METRICS.filter(function (m) { return m.col === 'median_depth' || m.col === 'total_reads'; })); } },
      { id: 'depth-hist', label: 'Target depth', build: buildDepthHistogramSubtab },
      { id: 'reference', label: 'Reference quality', build: function (r) { buildCohortMetricsSubtab(r, QC_METRICS.filter(function (m) { return m.col === 'max_corr' || m.col === 'n_refs'; })); } },
      { id: 'variation', label: 'Sample variability', build: function (r) { buildCohortMetricsSubtab(r, QC_METRICS.filter(function (m) { return m.col === 'noise_z' || m.col === 'n_cnvs'; })); } },
      { id: 'depth-corr', label: 'Depth vs. best correlation', build: buildDepthCorrelationSubtab },
      { id: 'pca', label: 'Cohort PCA', build: buildPcaSubtab }
    ], 'coverage');
    root.appendChild(st.el);
    onState(function () { st.refresh(); });
  });

  // Sample QC table is embedded in Overview > Samples.
  function buildSampleTable(root) {
    // Build a display table with both ECHO and CANOPE columns
    var dispCols = ['sample', 'status', 'level', 'n_cnvs', 'n_del', 'n_dup', 'n_high', 'n_medium', 'n_low',
                    'median_depth', 'total_reads', 'max_corr', 'n_refs', 'sex_inferred'];
    // Map ECHO column names
    var colMap = {};
    dispCols.forEach(function (c) {
      if (hasCol(S, c)) colMap[c] = c;
      else if (c === 'sample' && hasCol(S, 'Sample')) colMap[c] = 'Sample';
      else if (c === 'status' && hasCol(S, 'status')) colMap[c] = 'status';
      else if (c === 'level' && hasCol(S, 'Status')) colMap[c] = 'Status';
      else if (c === 'n_cnvs' && hasCol(S, 'CNV calls')) colMap[c] = 'CNV calls';
      else if (c === 'n_high' && hasCol(S, 'HIGH')) colMap[c] = 'HIGH';
      else if (c === 'n_medium' && hasCol(S, 'MEDIUM')) colMap[c] = 'MEDIUM';
      else if (c === 'n_low' && hasCol(S, 'LOW')) colMap[c] = 'LOW';
      else if (c === 'n_del' && hasCol(S, 'DEL')) colMap[c] = 'DEL';
      else if (c === 'n_dup' && hasCol(S, 'DUP')) colMap[c] = 'DUP';
      else if (c === 'median_depth' && hasCol(S, 'Median depth')) colMap[c] = 'Median depth';
      else if (c === 'total_reads' && hasCol(S, 'Total reads')) colMap[c] = 'Total reads';
      else if (c === 'max_corr' && hasCol(S, 'Max correlation')) colMap[c] = 'Max correlation';
      else if (c === 'n_refs' && hasCol(S, 'References')) colMap[c] = 'References';
    });
    var useCols = Object.keys(colMap);
    var T = { cols: useCols.map(function (c) { return colMap[c]; }), n: S.n, data: {} };
    useCols.forEach(function (c) { T.data[colMap[c]] = tcol(S, colMap[c]); });
    // add flags count column (CANOPE)
    if (IS_CANOPE) {
      var nf = SAMPLES.map(function (s) { var f = DATA.flags && DATA.flags[s]; return f ? f.filter(function (x) { return x.level === 'fail' || x.level === 'warn'; }).length : 0; });
      T.cols.push('flags'); T.data.flags = nf;
    }
    // Add sample-sheet attrs
    if (DATA.sample_attrs) Object.keys(DATA.sample_attrs).forEach(function (k) { T.cols.push(k); T.data[k] = asArray(DATA.sample_attrs[k]); });

    var tv = null;
    var labels = { n_cnvs: 'CNVs', n_del: 'DEL', n_dup: 'DUP', n_high: 'HIGH', n_medium: 'MEDIUM', n_low: 'LOW',
                   median_depth: 'median depth', total_reads: 'total reads', max_corr: 'best corr.', n_refs: 'refs', sex_inferred: 'sex (inferred)', flags: 'flags' };
    var cellClass = {};
    var render = {
      level: function (v) { return chip(v, LEVEL_LABEL[v] || v); },
      status: function (v) { return v === 'called' ? 'analysed' : v === 'Yes' ? 'analysed' : v === 'No' ? 'excluded' : v; },
      sample: function (v) { return h('strong', { text: v }); }
    };
    // numeric flag checks (CANOPE)
    if (IS_CANOPE) {
      var NF = {
        median_depth: function (v) { return isNum(v) && v < (TH.qc_min_cov || 0); },
        total_reads: function (v) { return isNum(v) && v < (TH.qc_min_total_reads || 0); },
        max_corr: function (v) { return isNum(v) && v < (TH.qc_min_corr || 0); },
        noise_z: function (v) { return isNum(v) && v > (TH.qc_zscore || 99); },
        n_refs: function (v) { return isNum(v) && v < 5; }
      };
      Object.keys(NF).forEach(function (k) { if (colMap[k]) cellClass[colMap[k]] = function (v) { return NF[k](v) ? 'warn' : ''; }; });
      cellClass.status = function (v) { return v === 'failed' || v === 'excluded' ? 'fail' : ''; };
    }
    var host = h('div');
    root.appendChild(card('Sample QC table', 'Select a row to filter the other report views; select it again or choose All samples to clear.', host));
    tv = tableView(host, T, { filename: 'samples', pageSize: 20, sort: ['flags', 'desc'], labels: labels, render: render, cellClass: cellClass,
      selected: function () { return state.sample ? SIDX[state.sample] : -1; },
      onRow: function (i) { setSample(SAMPLES[i] === state.sample ? '' : SAMPLES[i]); } });
    onState(function () { tv.redraw(); });
  }

  // ================================================================ CNV CALLS (with folded 4-panel viewer)
  function chrNum(c) {
    var s = String(c).replace(/^chr/i, '').toUpperCase(), n = parseInt(s, 10);
    if (s === 'X') return 23; if (s === 'Y') return 24; if (s === 'M' || s === 'MT') return 25;
    return isFinite(n) ? n : 99;
  }
  function callFlags(i) {
    var out = [];
    if (hasCol(CALLS, 'homology_gene') && tcol(CALLS, 'homology_gene')[i]) out.push('homology gene');
    if (hasCol(CALLS, 'bg_flag') && tcol(CALLS, 'bg_flag')[i]) out.push('background');
    if (!hasPlot(i)) out.push('no plot');
    return out.join(', ');
  }
  // Call-field helpers shared by the overview and call viewer
  function qsomeBfCell(i) {
    var q = tcol(CALLS, 'q_some')[i], bf = tcol(CALLS, 'bf')[i];
    if (IS_CANOPE) return orNa(q);
    if (IS_ECHO) return orNa(bf);
    return orNa(q !== undefined ? q : bf);
  }

  function callsTableDisplay() {
    // ECHO follows the reference report's explicit columns. CANOPE uses the
    // corresponding fields it actually provides; absent fields are omitted.
    var cols = [], data = {}, n = CALLS.n;
    function add(label, source) {
      if (!hasCol(CALLS, source)) return;
      cols.push(label);
      data[label] = tcol(CALLS, source);
    }
    add('ID', hasCol(CALLS, 'id') ? 'id' : 'ID');
    add('Sample', hasCol(CALLS, 'sample') ? 'sample' : 'Sample');
    add('Chr', hasCol(CALLS, 'chr') ? 'chr' : 'Chr');
    add('Start', hasCol(CALLS, 'start') ? 'start' : 'Start');
    add('End', hasCol(CALLS, 'end') ? 'end' : 'End');
    add('Gene', hasCol(CALLS, 'gene') ? 'gene' : 'Gene');
    add('Type', hasCol(CALLS, 'type') ? 'type' : 'Type');
    if (IS_ECHO) {
      add('Exons', 'exons');
      add('Nexons', 'nexons');
      add('Ratio', 'ratio');
      add('BF', 'bf');
      add('Confidence', 'confidence');
      add('Correlation', 'correlation');
      add('N.comp', 'n_comp');
      add('Calibration', 'calibration');
    } else {
      add('Targets', hasCol(CALLS, 'num_targ') ? 'num_targ' : 'nexons');
      add('MLCN', hasCol(CALLS, 'mlcn') ? 'mlcn' : 'cn_label');
      add('Q_SOME', 'q_some');
      add('Refs', hasCol(CALLS, 'num_refs') ? 'num_refs' : 'n_comp');
      add('Confidence', hasCol(CALLS, 'confidence') ? 'confidence' : 'Confidence');
      add('Calibration', hasCol(CALLS, 'calibration') ? 'calibration' : 'Calibration');
    }
    return { cols: cols, n: n, data: data };
  }

  function genomeCallHistogram(box) {
    var chrCol = hasCol(CALLS, 'chr') ? 'chr' : 'Chr';
    var chromosomes = uniq(tcol(CALLS, chrCol)).sort(function (a, b) { return chrNum(a) - chrNum(b); });
    if (!chromosomes.length) { box.appendChild(emptyNote('No chromosome coordinates were available for the calls.')); return; }
    plot(box, function (t) {
      var c = SEM[t], rows = range(CALLS.n).filter(function (i) { return !state.sample || callSample(i) === state.sample; });
      function counts(type) {
        return chromosomes.map(function (chr) {
          return rows.filter(function (i) { return String(tcol(CALLS, chrCol)[i]) === String(chr) && callType(i) === type; }).length;
        });
      }
      var del = counts('DEL'), dup = counts('DUP');
      var data = [
        { type: 'bar', name: 'DEL', x: chromosomes, y: del, marker: { color: c.DEL }, hovertemplate: '%{x}<br>DEL: %{y}<extra></extra>' },
        { type: 'bar', name: 'DUP', x: chromosomes, y: dup, marker: { color: c.DUP }, hovertemplate: '%{x}<br>DUP: %{y}<extra></extra>' }
      ];
      return { data: data, layout: baseLayout(t, { xlab: 'Chromosome', ylab: 'CNV calls', barmode: 'stack', legend: true,
        xaxis: { type: 'category', tickangle: -45 }, yaxis: { rangemode: 'tozero', dtick: 1 },
        margin: { l: 56, r: 16, t: 16, b: 76 } }) };
    });
  }

  // ---- 4-panel call viewer ----
  // ---- 4-panel call viewer ----
  function viewerTargetName(P, j) {
    var xlab = asArray(P.xlab), off = P.off ? asArray(P.off) : null;
    return off && off[j] ? xlab[j] : 'exon ' + (xlab[j] || '');
  }
  function viewerWindowTraces(P, t, y, refY, tip, refSize, showRefs) {
    var c = SEM[t], px = asArray(P.px), affected = asArray(P.affected), refNames = asArray(P.ref_names), data = [];
    if (showRefs && refY.length) {
      var rx = [], ry = [], rt = [];
      refY.forEach(function (vals, r) {
        vals.forEach(function (v, j) {
          rx.push(px[j]); ry.push(v); rt.push((refNames[r] || ('Reference ' + (r + 1))) + '<br>' + viewerTargetName(P, j));
        });
      });
      data.push({ type: 'scatter', mode: 'markers', name: 'Reference samples', x: rx, y: ry, text: rt,
        hoverinfo: 'text+y', marker: { color: c.ref, size: refSize, opacity: 0.7 } });
    }
    var pts = function (keep) {
      var ii = range(px.length).filter(function (j) { return !!affected[j] === keep; });
      return { x: ii.map(function (j) { return px[j]; }), y: ii.map(function (j) { return y[j]; }), text: ii.map(tip) };
    };
    var unaffected = pts(false), called = pts(true), hoverMode = showRefs && refY.length ? 'text+y' : 'text';
    data.push({ type: 'scatter', mode: 'markers', name: 'Test sample', x: unaffected.x, y: unaffected.y,
      text: unaffected.text, hoverinfo: hoverMode, marker: { color: c.test, size: 8 } });
    data.push({ type: 'scatter', mode: 'markers', name: 'Affected exon(s)', x: called.x, y: called.y,
      text: called.text, hoverinfo: hoverMode, marker: { color: c.hit, size: 11, symbol: 'diamond' } });
    return data;
  }
  function viewerPlots(P, cov, ci, z, showRefs) {
    var px = asArray(P.px), xlab = asArray(P.xlab), group = asArray(P.gene_group || P.group);
    var refCov = P.ref_cov ? asArray(P.ref_cov).map(asArray) : [];
    var refZ = P.ref_z ? asArray(P.ref_z).map(asArray) : [];
    var test = asArray(P.test_log || P.test), ratio = asArray(P.ratio), lo = asArray(P.lo), hi = asArray(P.hi), zScore = asArray(P.z);
    var axis = function () {
      var low = Math.min.apply(null, px) - 0.8, high = Math.max.apply(null, px) + 0.8;
      return { tickmode: 'array', tickvals: px, ticktext: xlab, tickangle: -45, range: [low, high], showgrid: false };
    };
    var tipCoverage = function (j) { return viewerTargetName(P, j); };
    plot(cov, function (t) {
      var data = viewerWindowTraces(P, t, test, refCov, tipCoverage, 6, showRefs);
      return { data: data, layout: baseLayout(t, { ylab: 'log(coverage)', legend: true, xaxis: axis() }) };
    });
    plot(ci, function (t) {
      var data = [], c = CHROME[t], groups = group.length ? uniq(group) : ['all'];
      groups.forEach(function (g, k) {
        var ii = range(px.length).filter(function (j) { return group.length ? group[j] === g : true; });
        if (!ii.length) return;
        var xs = ii.length === 1 ? [px[ii[0]] - 0.3, px[ii[0]] + 0.3] : ii.map(function (j) { return px[j]; });
        var band = function (values) { return ii.length === 1 ? [values[ii[0]], values[ii[0]]] : ii.map(function (j) { return values[j]; }); };
        data.push({ type: 'scatter', mode: 'lines', x: xs, y: band(lo), line: { width: 0 }, hoverinfo: 'skip', showlegend: false });
        data.push({ type: 'scatter', mode: 'lines', x: xs, y: band(hi), line: { width: 0 }, fill: 'tonexty', fillcolor: c.grid,
          name: '95% interval', showlegend: k === 0, hoverinfo: 'skip' });
      });
      var tipRatio = function (j) {
        return viewerTargetName(P, j) + '<br>observed ' + fmt(P.observed && P.observed[j]) +
          '<br>expected ' + fmt(P.expected && P.expected[j]) + '<br>ratio ' + fmt(ratio[j]) +
          ' [' + fmt(lo[j]) + ' – ' + fmt(hi[j]) + ']';
      };
      var points = viewerWindowTraces(P, t, ratio, [], tipRatio, 0, false);
      points[0].name = 'Observed'; points[0].marker.size = 9; points[1].name = 'Affected';
      data = data.concat(points);
      return { data: data, layout: baseLayout(t, { ylab: 'Observed / expected', legend: true, xaxis: axis(),
        shapes: [{ type: 'line', xref: 'paper', x0: 0, x1: 1, yref: 'y', y0: 1, y1: 1,
          line: { color: c.axis, width: 1, dash: 'dash' } }] }) };
    });
    plot(z, function (t) {
      var data = viewerWindowTraces(P, t, zScore, refZ, function (j) { return viewerTargetName(P, j); }, 5, showRefs);
      var c = CHROME[t];
      return { data: data, layout: baseLayout(t, { ylab: 'Z-score', legend: true, xaxis: axis(),
        shapes: [{ type: 'line', xref: 'paper', x0: 0, x1: 1, yref: 'y', y0: 0, y1: 0,
          line: { color: c.axis, width: 1, dash: 'dash' } }] }) };
    });
  }
  function renderGeneTrack(host, P) {
    clear(host);
    var px = asArray(P.px), genes = asArray(P.gene), groups = asArray(P.gene_group || P.group), off = asArray(P.off);
    if (!genes.length) {
      host.appendChild(emptyNote('Gene annotations are not available for this call.'));
      return;
    }
    var runs = [], begin = 0;
    function keyAt(i) {
      var grp = groups.length ? String(groups[i]) : '';
      return grp + '\u0000' + String(genes[i] || '');
    }
    for (var i = 1; i <= genes.length; i++) {
      if (i === genes.length || keyAt(i) !== keyAt(begin)) {
        var name = String(genes[begin] || '').trim() || 'Unannotated target';
        runs.push({ name: name, start: begin, end: i - 1, count: i - begin,
          off: !!off[begin], first: P.xlab && P.xlab[begin], last: P.xlab && P.xlab[i - 1] });
        begin = i;
      }
    }
    var track = h('div', { class: 'gene-track', role: 'group', 'aria-label': 'Genes across the selected CNV window' });
    runs.forEach(function (r, idx) {
      var cls = 'gene-track-segment track-' + ['a', 'b', 'c', 'd', 'e', 'f'][idx % 6] + (r.off ? ' off' : '');
      var range = r.first !== undefined && r.first !== null ? String(r.first) + (r.last !== undefined && r.last !== null && r.last !== r.first ? '–' + String(r.last) : '') : '';
      var meta = r.count + ' target' + (r.count === 1 ? '' : 's') + (range ? ' · ' + range : '');
      var segment = h('div', { class: cls, title: r.name + ' · ' + meta },
        h('div', { class: 'gene-track-name', text: r.name }),
        h('div', { class: 'gene-track-meta', text: meta }));
      segment.style.flexGrow = String(r.count);
      track.appendChild(segment);
    });
    host.appendChild(track);
  }


  tab('calls', 'CNV calls', function () { return true; }, function (root) {
    var stack = h('div', { class: 'stack' });
    if (!CALLS.n) {
      stack.appendChild(alertBox('info', 'No CNV was called.', 'A run without calls can be a valid result for the chromosomes analysed; see chromosome scope in Overview.'));
      root.appendChild(stack); return;
    }

    var confidenceValues = tcol(CALLS, hasCol(CALLS, 'confidence') ? 'confidence' : 'Confidence');
    var confLevels = ['HIGH', 'MEDIUM', 'LOW'];
    if (confidenceValues.some(function (v) { return v && confLevels.indexOf(String(v)) < 0; })) confLevels.push('NA');
    var f = { type: '', chr: '', conf: {} };
    confLevels.forEach(function (k) { f.conf[k] = true; });

    var typeSel = selectOf([
      { value: '', text: 'All types' },
      { value: 'DEL', text: 'Deletions' },
      { value: 'DUP', text: 'Duplications' }
    ], '', null, 'Call type');
    var chrCol = hasCol(CALLS, 'chr') ? 'chr' : 'Chr';
    var chrs = uniq(tcol(CALLS, chrCol)).sort(function (a, b) { return chrNum(a) - chrNum(b); });
    var chrSel = selectOf([{ value: '', text: 'All chromosomes' }].concat(chrs.map(function (c) {
      return { value: c, text: c };
    })), '', null, 'Chromosome');

    var confBox = h('div', { class: 'seg', role: 'group', 'aria-label': 'Confidence' });
    confLevels.forEach(function (level) {
      var label = level === 'NA' ? 'Unscored' : level;
      confBox.appendChild(h('button', {
        type: 'button', 'aria-pressed': 'true', text: label,
        onclick: function () {
          f.conf[level] = !f.conf[level];
          this.setAttribute('aria-pressed', f.conf[level] ? 'true' : 'false');
          apply();
        }
      }));
    });

    var count = h('span', { class: 'muted small' });
    var tableHost = h('div'), detail = h('div', { class: 'call-detail-wrap' });
    var Tdisp = callsTableDisplay();
    var table = null;

    function pass(i) {
      if (state.sample && callSample(i) !== state.sample) return false;
      if (f.type && callType(i) !== f.type) return false;
      if (f.chr && String(tcol(CALLS, chrCol)[i]) !== f.chr) return false;
      var ck = callConf(i);
      ck = (ck === 'HIGH' || ck === 'MEDIUM' || ck === 'LOW') ? ck : 'NA';
      return f.conf[ck] !== false;
    }

    stack.appendChild(h('div', { class: 'toolbar calls-toolbar' },
      field('Confidence', confBox), field('Type', typeSel), field('Chromosome', chrSel), count));
    stack.appendChild(card('CNV calls', 'Select a row to inspect its coverage evidence below.', tableHost));
    stack.appendChild(detail);
    root.appendChild(stack);

    table = tableView(tableHost, Tdisp, {
      filename: 'echo_cnv_calls', pageSize: 20,
      sort: ['ID', 'asc'],
      wrap: ['Gene'],
      render: {
        Type: typeChip,
        Confidence: confChip,
        Sample: function (v) { return h('strong', { text: v }); },
        Calibration: function (v) { return v === 'check' ? chip('warn', 'check') : v ? chip('ok', v) : ''; },
        Start: function (v) { return fmtInt(v); },
        End: function (v) { return fmtInt(v); }
      },
      selected: function () { return state.call === null || state.call === '' ? -1 : CALL_ROW[state.call]; },
      onRow: function (i) {
        var clicked = callId(i);
        state.call = state.call === clicked ? '' : clicked;
        notify();
      }
    });

    function apply() {
      table.setPredicate(pass);
      var visible = table.visible();
      count.textContent = visible.length + ' of ' + CALLS.n + ' calls';
      if (state.call !== null && state.call !== '' && visible.indexOf(CALL_ROW[state.call]) < 0) state.call = null;
      if (state.call === null && visible.length) state.call = callId(visible[0]);
      var row = (state.call === null || state.call === '') ? undefined : CALL_ROW[state.call];
      if (row !== undefined && visible.indexOf(row) >= 0) showCallDetail(detail, row);
      else {
        clear(detail);
        detail.appendChild(card('Call evidence', null, emptyNote('No call selected. Select a row to display its values and figures; select the same row again to clear the selection.')));
      }
      // Redraw once the current call is known so the selected table row is visible.
      table.redraw();
    }

    typeSel.addEventListener('change', function () { f.type = typeSel.value; apply(); });
    chrSel.addEventListener('change', function () { f.chr = chrSel.value; apply(); });
    onState(apply);
    apply();
  });

  function showCallDetail(host, i) {
    clear(host);
    var id = callId(i), P = PLOTS[String(id)];
    var head = h('div', { class: 'detail-head' },
      h('h2', { text: callSample(i) + ' \u00b7 ' + callGene(i) }),
      typeChip(callType(i)), confChip(callConf(i)),
      h('span', { class: 'muted', text: callInterval(i) }));
    // Only add model-level values that are not already shown in the calls table.
    var kvPairs = [], additionalFlags = callFlags(i);
    if (additionalFlags) kvPairs.push(['Additional flags', additionalFlags]);
    if (IS_ECHO) {
      kvPairs = kvPairs.concat([
        ['Observed count', orNa(tcol(CALLS, 'observed')[i])],
        ['Expected count', orNa(tcol(CALLS, 'expected')[i])],
        ['Comparators', orNa(tcol(CALLS, 'comparators')[i])]
      ]);
    }
    if (IS_CANOPE) {
      kvPairs = kvPairs.concat([
        ['KB', orNa(tcol(CALLS, 'kb')[i])],
        ['Transition / interval details', orNa(tcol(CALLS, 'transition_details')[i])]
      ]);
    }
    var body = h('div', { class: 'stack' });
    if (!P) {
      body.appendChild(emptyNote(META.windows_limited ? 'Coverage profiles were only embedded for the first ' + META.windows_limit + ' calls.' : 'No coverage window for this call.'));
    } else {
      // Call-level calibration is shown in the selected call details and call table.
      var stripBox = h('div', { class: 'gene-track-container' }), covBox = plotBox('mid'), ciBox = plotBox('mid'), zBox = plotBox('mid');
      var refsChk = h('input', { type: 'checkbox', checked: true, 'aria-label': 'Show reference samples' });
      body.appendChild(h('div', { class: 'toolbar' }, h('label', { class: 'field inline' }, refsChk, 'Show reference samples')));
      body.appendChild(stripBox);
      body.appendChild(card('Log coverage (scaled references vs test sample)', null, covBox));
      body.appendChild(card(IS_ECHO ? 'Observed / expected with 95% beta-binomial interval' : 'Observed / expected with 95% negative-binomial interval', null, ciBox));
      body.appendChild(card('Z-score vs references', null, zBox));
      function drawPanels() {
        renderGeneTrack(stripBox, P);
        viewerPlots(P, covBox, ciBox, zBox, refsChk.checked);
      }
      refsChk.addEventListener('change', drawPanels);
      drawPanels();
    }
    var detailGrid = kvPairs.length ? h('div', { class: 'grid g-main' }, h('div', { class: 'stack' }, body), kv(kvPairs)) : body;
    host.appendChild(card('Call details', 'CNV.ID ' + id, h('div', null, head, detailGrid)));
    applyGlobalAlertFilter();
    host.scrollIntoView && host.scrollIntoView({ block: 'nearest', behavior: 'smooth' });
  }

  function pcaPlot(div) {
    var P = PCA;
    plot(div, function (t) {
      var names = P.names || P.samples, pc1 = P.pc1 || P.pcs.PC1, pc2 = P.pc2 || P.pcs.PC2;
      var groups = P.group ? uniq(P.group) : ['All'], cc = catColors(Math.max(groups.length, 1), t), s = state.sample;
      var data = groups.map(function (g, k) {
        var ii = range(names.length).filter(function (i) { return (P.group ? P.group[i] : 'All') === g; });
        return { type: 'scatter', mode: 'markers+text', name: g, x: ii.map(function (i) { return pc1[i]; }), y: ii.map(function (i) { return pc2[i]; }),
          text: ii.map(function (i) { return names[i]; }), textposition: 'top right', textfont: { size: 10, color: CHROME[t].secondary },
          hovertemplate: '<b>%{text}</b><br>PC1 %{x:.2f}<br>PC2 %{y:.2f}<extra>' + g + '</extra>',
          marker: { color: cc[k], size: ii.map(function (i) { return s && names[i] === s ? 15 : 9; }), line: { width: ii.map(function (i) { return s && names[i] === s ? 3 : 0; }), color: CHROME[t].text } } };
      });
      var v1 = P.var1 || (P.var ? P.var[0] : 0), v2 = P.var2 || (P.var ? P.var[1] : 0);
      return { data: data, layout: baseLayout(t, { xlab: 'PC1 (' + fmt(v1) + '%)', ylab: 'PC2 (' + fmt(v2) + '%)', legend: groups.length > 1, margin: { l: 56, r: 24, t: 28, b: 48 } }) };
    });
  }
  function depthHist(div) {
    var Hh = DEPTH_HIST;
    plot(div, function (t) {
      var sem = SEM[t];
      var data = [{ type: 'bar', x: Hh.mid, y: Hh.counts, marker: { color: Hh.mid.map(function (m) { return (TH.qc_min_cov || TH.min_cov) && Math.pow(10, m) < (TH.qc_min_cov || TH.min_cov) ? sem.bad : sem.dup; }) },
        hovertemplate: 'median depth \u2248 %{customdata:.0f}<br>%{y} targets<extra></extra>', customdata: Hh.mid.map(function (m) { return Math.pow(10, m); }) }];
      var thr = Math.log10(Math.max(TH.qc_min_cov || TH.min_cov || 1, 1));
      return { data: data, layout: baseLayout(t, { xlab: 'Median depth per target (log)', ylab: 'Targets',
        shapes: (TH.qc_min_cov || TH.min_cov) ? [{ type: 'line', xref: 'x', yref: 'paper', x0: thr, x1: thr, y0: 0, y1: 1, line: { color: sem.warn, width: 2, dash: 'dash' } }] : [] }) };
    });
  }

  function buildChromosomeOverview(root) {
    var scopeCard = chromosomeScopeCard();
    if (scopeCard) { scopeCard.id = 'chromosome-scope'; root.appendChild(scopeCard); }
    var inv = CHR.inventory;
    if (inv && inv.n) {
      var invHost = h('div'), barBox = plotBox('short');
      tableView(invHost, inv, { filename: 'chromosomes', pageSize: 30,
        labels: { targets_in_bed: 'targets in BED', targets_analysed: 'targets analysed' },
        render: { analysed: function (v) { return v ? chip('ok', 'analysed') : chip('warn', 'not analysed'); } },
        cellClass: { targets_analysed: function (v, i) { return tcol(inv, 'analysed')[i] ? '' : 'warn'; } } });
      root.appendChild(card('Chromosome inventory', 'Targets in BED and analysed status by chromosome', invHost));
      root.appendChild(card('Targets per chromosome', 'Targets in BED vs. analysed', barBox));
      plot(barBox, function (t) {
        var c = tcol(inv, 'chromosome');
        return { data: [
          { type: 'bar', name: 'in BED', x: c, y: tcol(inv, 'targets_in_bed'), marker: { color: SEM[t].NA } },
          { type: 'bar', name: 'analysed', x: c, y: tcol(inv, 'targets_analysed'), marker: { color: CAT[t][0] } }
        ], layout: baseLayout(t, { ylab: 'targets', barmode: 'overlay', legend: true,
          xaxis: { type: 'category', tickangle: -45 }, margin: { l: 56, r: 16, t: 10, b: 80 } }) };
      });
    }
    var D = CHR.depth;
    if (D && D.chromosomes && asArray(D.chromosomes).length) {
      var hbox = plotBox('tall');
      root.appendChild(card('Chromosome depth relative to autosomes', 'Median target depth / sample autosomal median', hbox));
      plot(hbox, function (t) {
        var dz = asArray(D.z).map(asArray);
        return { data: [{ type: 'heatmap', x: asArray(D.samples), y: asArray(D.chromosomes), z: dz, zmin: 0, zmax: 2,
          colorscale: [[0, DIV[0]], [0.25, DIV[1]], [0.5, DIV[2]], [0.75, DIV[3]], [1, DIV[4]]],
          hovertemplate: '%{x}<br>%{y}: %{z:.2f}<extra></extra>', colorbar: { title: { text: 'ratio' } } }],
          layout: baseLayout(t, { xaxis: { type: 'category', tickangle: -45 }, yaxis: { type: 'category', autorange: 'reversed' }, margin: { l: 70, r: 16, t: 10, b: 90 } }),
          onClick: function (ev) { var p = ev.points && ev.points[0]; if (p) setSample(p.x); } };
      });
    }
    var X = CHR.sex;
    if (X && X.n) {
      var sbox = plotBox(''), sh = h('div');
      tableView(sh, X, { filename: 'sex_check', pageSize: 10,
        labels: { x_ratio: 'chrX ratio', y_ratio: 'chrY ratio' },
        cellClass: { mismatch: function (v) { return v ? 'warn' : ''; } },
        onRow: function (i) { setSample(tcol(X, 'sample')[i]); } });
      root.appendChild(h('div', { class: 'grid g2' },
        card('Sex check', 'chrX vs chrY depth ratio', sbox), card('Sex check table', null, sh)));
      plot(sbox, function (t) {
        var names = asArray(X.sample), xr = asArray(X.x_ratio), yr = asArray(X.y_ratio);
        return { data: [{ type: 'scatter', mode: 'markers+text', x: xr, y: yr, text: names, textposition: 'top right',
          marker: { color: asArray(X.mismatch).map(function (v) { return v ? SEM[t].warn : SEM[t].ok; }), size: 9 },
          hovertemplate: '%{text}<br>chrX ratio %{x}<br>chrY ratio %{y}<extra></extra>' }],
          layout: baseLayout(t, { xlab: 'chrX ratio', ylab: 'chrY ratio' }),
          onClick: function (ev) { var p = ev.points && ev.points[0]; if (p && p.text) setSample(p.text); } };
      });
    }
  }

  // ================================================================ METHODS
  tab('methods', 'Methods', function () { return true; }, function (root) {
    var stack = h('div', { class: 'stack' });
    var R = CONF_RULES, sc = SCORING;

    // 1. Confidence score (both ECHO scoring and CANOPE conf_rules)
    if (sc || R) {
      var confRulesList = h('ul', null,
        h('li', null, h('b', { text: 'HIGH' }), ' only if all HIGH criteria are met.'),
        h('li', null, h('b', { text: 'MEDIUM' }), ' if HIGH is not met but all MEDIUM criteria are.'),
        h('li', null, h('b', { text: 'LOW' }), ' otherwise.'));
      var lowGenes = [], seenLow = {};
      [].concat(R && R.genes ? asArray(R.genes) : [], sc && sc.forced_low ? asArray(sc.forced_low) : []).forEach(function (g) {
        var gene = String(g || '').trim(), key = gene.toUpperCase();
        if (gene && !seenLow[key]) { seenLow[key] = true; lowGenes.push(gene); }
      });
      if (lowGenes.length) confRulesList.appendChild(h('li', null, h('b', { text: 'Always LOW' }), ' (homology): ' + lowGenes.join(', ')));
      if (sc && sc.terminal) confRulesList.appendChild(h('li', { text: sc.terminal }));
      var confTable;
      if (sc && sc.rows) {
        confTable = h('table', { class: 'tbl' },
          h('thead', null, h('tr', null, ['Category', 'Correlation', 'References', 'Fold change'].map(function (c) { return h('th', { text: c }); }))),
          h('tbody', null, sc.rows.map(function (r) { return h('tr', null, r.map(function (v, j) { return h('td', null, j === 0 ? chip(v, v) : v); })); })));
      } else {
        confTable = h('table', { class: 'tbl' },
          h('thead', null, h('tr', null, ['Category', 'Q_SOME / BF', 'References', 'Targets'].map(function (c) { return h('th', { text: c }); }))),
          h('tbody', null,
            h('tr', null, h('td', null, chip('conf HIGH', 'HIGH')), h('td', { text: '>= ' + (R.high_q || '') }), h('td', { text: '>= ' + (R.high_refs || '') }), h('td', { text: '>= ' + (R.high_targ || '') })),
            h('tr', null, h('td', null, chip('conf MEDIUM', 'MEDIUM')), h('td', { text: '>= ' + (R.med_q || '') }), h('td', { text: '>= ' + (R.med_refs || '') }), h('td', { text: '>= ' + (R.med_targ || '') })),
            h('tr', null, h('td', null, chip('conf LOW', 'LOW')), h('td', { text: 'anything else' }), h('td', { text: 'any' }), h('td', { text: 'any' }))));
      }
      stack.appendChild(detailCard('Confidence score', null, h('div', null, confRulesList, h('div', { class: 'tbl-wrap', style: { marginTop: '10px' } }, confTable))));
    }

    // 2. Reading the plots
    stack.appendChild(detailCard('Reading the plots', null, h('div', { class: 'prose' },
      h('p', { text: 'Each CNV call shows the test sample against its selected reference samples over the affected exons ± 5 neighbouring targets. Coverage is shown on a log scale after scaling each reference to the median coverage of the test sample in the window.' }),
      h('p', { text: IS_ECHO
        ? 'Observed / expected is the test count divided by the count expected from the reference aggregate; the grey ribbon is the 95% interval of the beta-binomial model (overdispersion phi fitted by ExomeDepth). The z-score uses the model’s own standard deviation, so reference samples should scatter around 0.'
        : 'Observed / expected is the test count divided by the count expected from the reference aggregate; the grey ribbon is the 95% predictive interval of the negative-binomial model. The z-score is measured against the reference distribution, so reference samples should scatter around 0.' }),
      h('p', { text: 'If more background exons than expected fall outside the interval, the call is flagged “check” (poorly calibrated reference set).' }))));
    // 4. Calling model
    stack.appendChild(detailCard(IS_ECHO ? 'How ECHO works' : 'How CANOPE works', null, h('div', { class: 'prose' },
      IS_ECHO
        ? [h('p', { text: 'ECHO is an ExomeDepth-based exome CNV caller. For every sample it selects the most correlated reference samples, combines them into an expected coverage per exon, and models the read count with a beta-binomial distribution whose overdispersion (phi) is fitted from the references. A Bayes factor scores each candidate call; confidence (HIGH / MEDIUM / LOW) is derived from the correlation, the number of references and the fold change.' }),
           h('p', { text: 'One run covers one chromosome mode: autosomes (A), chrX (X) or chrY (Y). Sex chromosomes are therefore never analysed together with the autosomes.' })]
        : [h('p', { text: 'CANOPE is a CANOES-derived exome CNV caller. For every sample it selects the most correlated reference samples, combines them with non-negative weights into an expected coverage per exon, and models the read count with a negative binomial whose variance is estimated from the references. A hidden Markov model with distance-decayed transitions (Viterbi) calls deletions and duplications; Q_SOME is the phred-scaled probability that at least one exon of the call is abnormal, and MLCN the most likely copy number.' }),
           h('p', { text: 'One run covers one chromosome mode: autosomes (A), chrX with female references (XX) or chrX/chrY with male references (XY).' })])));

    // 5. Warning glossary
    var glossary = [
      ['Chromosomes not analysed', 'The BED/count table contains targets that this run did not call. A missing call there is not a negative result.'],
      ['Sample failed / excluded', 'Failed samples raised an error while calling; excluded samples were dropped before calling. Neither has CNV calls, and excluded samples are not used as references.'],
      ['Low depth / low reads', 'Median target depth or total reads below the QC threshold: reduced sensitivity, especially for single-exon events.'],
      ['Low correlation / few references', 'No other sample resembles this one (batch, capture kit, tissue), so the expected coverage is poor and false calls are likely.'],
      ['Noisy sample', 'The exon-count dispersion is an outlier within the cohort (robust z-score).'],
      ['Genotyping skipped', 'Too many calls were produced; Q_SOME and MLCN were not computed and every call is LOW.'],
      ['Sex / mode mismatch', 'The annotated sex disagrees with chrX/chrY coverage, or the sample does not match the sex expected by the mode.'],
      ['Background calibration', 'A high share of non-called exons around a call fall outside their 95% interval, meaning the expected-coverage model is miscalibrated.']];
    stack.appendChild(detailCard('Warning glossary', null, kv(glossary)));

    root.appendChild(stack);
  });

  // ================================================================ RUN
  tab('run', 'Run', function () { return true; }, function (root) {
    var L = LOG || {}, blocks = [];
    var summary = runKv();
    if (summary.childElementCount) blocks.push(detailCard('Report metadata', null, summary));

    // Show the original configuration OR a resolved fallback, never both.
    if (DATA.config_yaml) blocks.push(detailCard('Configuration used', 'Original configuration for this run', h('pre', { text: DATA.config_yaml })));
    else if (SETTINGS && SETTINGS.length) blocks.push(detailCard('Resolved settings', null, kv(SETTINGS.map(function (s) { return [s[0], s[1]]; }))));
    else if (Object.keys(PARAMS).length) blocks.push(detailCard('Run parameters', null,
      kv(Object.keys(PARAMS).map(function (k) { return [k, Array.isArray(PARAMS[k]) ? PARAMS[k].join(', ') : String(PARAMS[k])]; }))));

    if (DATA.sample_table && DATA.sample_table.n) {
      var sampleManifest = h('div');
      tableView(sampleManifest, DATA.sample_table, { filename: 'sample_manifest', pageSize: 15 });
      blocks.push(detailCard('Input sample manifest', null, sampleManifest));
    }
    if (DATA.ref_bams && DATA.ref_bams.n) {
      var refs = h('div');
      tableView(refs, DATA.ref_bams, { filename: 'reference_bams', pageSize: 15 });
      blocks.push(detailCard('External reference BAMs', null, refs));
    }
    // Keep the raw execution log in one place; structured QC flags remain in Samples.
    if (L.tail && L.tail.length) blocks.push(detailCard('Pipeline log' + (L.truncated ? ' (last lines)' : ''), null, h('pre', { text: L.tail.join('\n') })));
    if (L.session && L.session.length) blocks.push(detailCard('Session information', null, h('pre', { text: L.session.join('\n') })));

    if (!blocks.length) {
      root.appendChild(emptyNote('No run metadata, configuration or execution log was recorded.'));
      return;
    }
    root.appendChild(h('div', { class: 'stack' }, blocks));
  });

  // ================================================================ shell
  var panels = {}, tabButtons = {}, active = null;
  function activate(id, push) {
    var t = TABS.filter(function (x) { return x.id === id && x.available(); })[0] || TABS.filter(function (x) { return x.available(); })[0];
    if (!t) return;
    id = t.id;
    Object.keys(panels).forEach(function (k) {
      var on = k === id; panels[k].hidden = !on;
      tabButtons[k].setAttribute('aria-selected', on ? 'true' : 'false'); tabButtons[k].tabIndex = on ? 0 : -1;
    });
    if (!t.built) {
      t.built = true;
      try { t.build(panels[id]); }
      catch (e) { panels[id].appendChild(h('div', { class: 'empty', text: 'This tab could not be built: ' + e.message })); if (window.console) console.error(e); }
    } else refreshVisible(panels[id]);
    applyGlobalAlertFilter();
    setTimeout(function () { refreshVisible(panels[id]); applyGlobalAlertFilter(); }, 60);
    active = id;
    if (push) { try { history.replaceState(null, '', '#tab=' + id); } catch (e) {} }
  }
  function init() {
    applyThemeAttr();
    var reportTitle = String(META.title || 'ECHO report').replace(/^ECHO\s+report\s*-\s*(.+)$/i, 'ECHO - $1 report');
    document.title = reportTitle;
    var tabsEl = $('#tabs'), main = $('#panels');
    var tabOrder = { overview: 0, samples: 1, qc: 2, calls: 3, methods: 4, run: 5 };
    var avail = TABS.filter(function (t) { return t.available(); }).sort(function (a, b) {
      var ao = tabOrder[a.id] === undefined ? 99 : tabOrder[a.id];
      var bo = tabOrder[b.id] === undefined ? 99 : tabOrder[b.id];
      return ao - bo;
    });
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
    // header sample picker
    var pick = $('#sample-pick');
    if (pick && SAMPLES.length > 0) {
      headerPicker = h('select', { id: 'sample-sel', 'aria-label': 'Filter report by sample' },
        h('option', { value: '', text: 'All samples (' + SAMPLES.length + ')' }),
        SAMPLES.map(function (s) { return h('option', { value: s, text: s }); }));
      headerPicker.value = state.sample;
      headerPicker.addEventListener('change', function () { setSample(headerPicker.value); });
      pick.appendChild(h('label', { for: 'sample-sel', text: 'Filter sample' }));
      pick.appendChild(headerPicker);
    } else if (pick) { pick.hidden = false; }
    if (pick) {
      globalSeverityControl = makeSeveritySegments(state.alertFilter, [
        { value: 'all', label: 'All' }, { value: 'fail', label: 'Fail' },
        { value: 'warning', label: 'Warn' }, { value: 'note', label: 'Note' }, { value: 'hide', label: 'Off' }
      ], setAlertFilter, 'Global alert severity filter');
      pick.appendChild(h('div', { class: 'global-alert-filter', title: 'Show or hide warning levels throughout the report' },
        h('span', { text: 'Severity' }), globalSeverityControl.el));
    }

    // Persistent chromosome-scope warning: this lives in the sticky shell so it
    // remains visible above every tab and is controlled by the single Severity filter.
    var noticeHost = $('#global-notice');
    var missingClasses = missingBedChromosomeClasses();
    if (noticeHost && missingClasses.length) {
      var scopeLink = h('button', { class: 'btn small-btn', type: 'button', text: 'Chromosome scope details', onclick: function () {
        activate('overview', true);
        if (overviewSubtabs) overviewSubtabs.switchTo('chromosomes');
      } });
      noticeHost.appendChild(alertBox('fail', 'Chromosome class missing from the BED file', missingClasses.join(', ') + ' not present in BED file.', scopeLink));
    }
    if (noticeHost) {
      var sampleData = S.data || {}, sampleNames = sampleData.sample || sampleData.Sample || SAMPLES;
      var analysed = sampleData.Analysed || sampleData.analysed || [];
      var statuses = sampleData.status || [];
      var notAnalysedNames = sampleNames.filter(function (name, i) {
        return (analysed.length && /^(no|false|0)$/i.test(String(analysed[i]))) ||
          (!analysed.length && /^(excluded|failed|not analysed|not analyzed)$/i.test(String(statuses[i] || '')));
      });
      if (notAnalysedNames.length) {
        noticeHost.appendChild(alertBox('fail', 'Samples not analysed (' + notAnalysedNames.length + ')',
          notAnalysedNames.join(', ') + ' not analysed.'));
      }
    }

    var btn = $('#theme-btn');
    var labelFn = function () { btn.textContent = 'Theme: ' + themePref; };
    labelFn();
    btn.addEventListener('click', function () {
      themePref = themePref === 'auto' ? 'light' : themePref === 'light' ? 'dark' : 'auto';
      store('unified-report-theme', themePref); applyThemeAttr(); labelFn(); rerenderAll();
    });
    if (window.matchMedia) {
      var mq = window.matchMedia('(prefers-color-scheme: dark)');
      var onChange = function () { if (themePref === 'auto') rerenderAll(); };
      if (mq.addEventListener) mq.addEventListener('change', onChange);
    }
    var m = /tab=([\w-]+)/.exec(location.hash);
    activate(m ? m[1] : 'overview', false);
    applyGlobalAlertFilter();
    window.addEventListener('hashchange', function () { var mm = /tab=([\w-]+)/.exec(location.hash); if (mm && mm[1] !== active) activate(mm[1], false); });
    document.documentElement.setAttribute('data-ready', '1');
  }

  function ensurePlotly(cb, timeoutMs) {
    timeoutMs = timeoutMs || 8000;
    var done = false;
    function has() { return (typeof window.Plotly !== 'undefined') && window.Plotly && (typeof window.Plotly.react === 'function'); }
    function finish(ok) {
      if (done) return; done = true; clearTimeout(timer);
      if (ok && has()) cb();
      else {
        var host = document.getElementById('panels');
        if (host) host.innerHTML = '<div class="empty">Plotly.js failed to load. Check your connection or reinstall the R package.</div>';
        if (window.console) console.error('[report] Plotly unavailable.');
      }
    }
    var timer = setTimeout(function () { finish(false); }, timeoutMs);
    if (has()) { finish(true); return; }
    var s = document.createElement('script');
    s.src = 'https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js';
    s.async = false; s.onload = function () { finish(true); }; s.onerror = function () { finish(false); };
    document.head.appendChild(s);
  }

  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', function () { ensurePlotly(init); });
  else ensurePlotly(init);
})();
)---"


# ----------------------------------------------------------------------------
# R helpers (shared by both packages)
# ----------------------------------------------------------------------------

#' CSS string for the unified report
#' @export
report_css <- function() REPORT_CSS

#' JS string for the unified report
#' @export
report_js <- function() REPORT_JS

#' Column-oriented table serialiser (data.frame -> list for JSON)
#' @export
rep_table <- function(df) {
  if (is.null(df)) return(list(cols = I(character()), n = 0L, data = list()))
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  data <- lapply(df, function(v) {
    if (is.factor(v)) v <- as.character(v)
    if (is.double(v)) v <- round(v, 4)
    I(unname(v))
  })
  list(cols = I(colnames(df)), n = nrow(df), data = data)
}

#' Empty table placeholder
#' @export
rep_empty_table <- function(cols) {
  rep_table(stats::setNames(rep(list(character()), length(cols)), cols))
}

#' Inert JSON for a <script type="application/json"> element
#' @export
report_json <- function(x) {
  json <- as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", na = "null",
                                        digits = NA, force = TRUE))
  json <- gsub("</", "<\\/", json, fixed = TRUE)
  gsub("<!--", "<\\!--", json, fixed = TRUE)
}

#' Read the plotly.js bundle from the installed plotly R package (or NULL)
#' @export
report_plotly_js <- function() {
  if (!requireNamespace("plotly", quietly = TRUE)) return(NULL)
  deps <- tryCatch(suppressWarnings(suppressMessages(htmltools::findDependencies(plotly::plot_ly()))),
                   error = function(e) list())
  dep <- Filter(function(d) identical(d$name, "plotly-main"), deps)
  if (!length(dep)) return(NULL)
  dep <- dep[[1]]
  path <- file.path(dep$src$file, dep$script)
  if (!file.exists(path)) path <- system.file(dep$src$file, dep$script, package = dep$package %||% "plotly")
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  raw <- readBin(path, "raw", file.size(path))
  txt <- rawToChar(raw); Encoding(txt) <- "UTF-8"
  gsub("</script", "<\\/script", txt, fixed = TRUE)
}

#' HTML-escape a string
#' @export
report_esc <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

#' Assemble the full HTML document from a payload and write it to out_file.
#'
#' Both echo_report.R and canope_report.R call this with their payload
#' (a list with a \code{meta$source} field set to "ECHO" or "CANOPE").
#' The resulting HTML is byte-identical in structure and styling regardless
#' of source; only the embedded JSON blob differs.
#'
#' @param payload List (the report data, with meta$source set).
#' @param settings Optional list; \code{settings$report_plotly} can be "inline" (default) or "cdn".
#' @param out_file Path to write the HTML file to.
#' @return Invisibly, out_file.
#' @export
report_html <- function(payload, settings = list(), out_file) {
  css <- report_css()
  js  <- report_js()
  mode <- tolower(as.character(settings$report_plotly %||% "inline")[1])
  plotly_tag <- NULL
  if (mode != "cdn") {
    src <- report_plotly_js()
    if (is.null(src)) {
      message("[WARNING] plotly.js bundle not found; loading from CDN (needs internet).")
    } else {
      plotly_tag <- paste0("<script>\n", src, "\n</script>\n")
    }
  }
  if (is.null(plotly_tag))
    plotly_tag <- "<script src=\"https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js\"></script>\n"

  title_value <- as.character(payload$meta$title %||% 'ECHO report')
  title_value <- sub('^ECHO\\s+report\\s*-\\s*(.+)$', 'ECHO - \\1 report', title_value, ignore.case = TRUE)
  payload$meta$title <- title_value
  title <- report_esc(title_value)
  html <- paste0(
    "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n",
    "<title>", title, "</title>\n<style>\n", css, "\n</style>\n</head>\n<body>\n",
    "<header class=\"top\"><div class=\"top-inner\">\n",
    "<div class=\"title-row\"><h1>", title, "</h1>",
    "<span class=\"spacer\"></span>",
    "<div id=\"sample-pick\" class=\"sample-pick\"></div>",
    "<button id=\"theme-btn\" class=\"btn\" type=\"button\" aria-label=\"Colour theme\">Theme</button></div>\n",
    "<div id=\"global-notice\" class=\"global-notice\" aria-live=\"polite\"></div>\n",
    "<nav id=\"tabs\" class=\"tabs\" role=\"tablist\" aria-label=\"Report sections\"></nav>\n",
    "</div></header>\n<main id=\"panels\"></main>\n",
    "<noscript><p style=\"padding:16px\">This report needs JavaScript enabled.</p></noscript>\n",
    plotly_tag,
    "<script id=\"unified-data\" type=\"application/json\">", report_json(payload), "</script>\n",
    "<script>\n", js, "\n</script>\n</body>\n</html>\n")
  con <- file(out_file, "wb")
  on.exit(close(con), add = TRUE)
  writeBin(charToRaw(enc2utf8(html)), con)
  message("[INFO] ", payload$meta$source, " report written: ", out_file, " (",
          format(structure(file.size(out_file), class = "object_size"), units = "auto"), ")")
  invisible(out_file)
}

# null-coalesce helper
`%||%` <- function(a, b) if (is.null(a)) b else a
