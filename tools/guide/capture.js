// Capture sections of the pipeline's run report as PNG figures for the user guide.
// Usage: node tools/guide/capture.js <run_report.html> <figure folder>
// Each shot clips from a heading to the next heading of the same or higher level (or to a named one).
const { chromium } = require('playwright');
const path = require('path');

const SHOTS = [
  ['Demo Road', 'rep_summary', 900, 'Plan'],
  ['Plan', 'rep_plan', 1500, 'Crossing schedule'],
  ['Crossing schedule', 'rep_schedule', 1100],
  ['DEM checks', 'rep_dem_checks', 1000, 'DEM uncertainty at crossings'],
  ['DEM uncertainty at crossings', 'rep_uncertainty', 1000, 'Land-cover scenario'],
  ['Land-cover scenario', 'rep_scenario', 800],
  ['Check against mapped drainage', 'rep_mapped', 1000, 'Time of concentration'],
  ['Time of concentration', 'rep_tc', 1100, 'Drainage coverage along the road'],
  ['Drainage coverage along the road', 'rep_coverage', 900],
  ['Inputs and provenance', 'rep_provenance', 1100],
  ['Warnings', 'rep_warnings', 600],
];

(async () => {
  const [report, out] = process.argv.slice(2).map((p) => path.resolve(p));
  const browser = await chromium.launch({ executablePath: '/opt/pw-browsers/chromium' }).catch(() => chromium.launch());
  const page = await browser.newPage({ viewport: { width: 1000, height: 1400 }, deviceScaleFactor: 2 });
  await page.goto('file://' + report, { waitUntil: 'networkidle' });
  for (const [heading, name, maxHeight, until] of SHOTS) {
    const box = await page.evaluate(([text, limit, end]) => {
      const heads = [...document.querySelectorAll('h1,h2,h3')];
      const index = heads.findIndex((h) => h.textContent.trim().startsWith(text));
      if (index < 0) return null;
      const head = heads[index];
      const level = Number(head.tagName[1]);
      const next = end
        ? heads.slice(index + 1).find((h) => h.textContent.trim().startsWith(end))
        : heads.slice(index + 1).find((h) => Number(h.tagName[1]) <= level);
      const top = head.getBoundingClientRect().top + window.scrollY - 8;
      const bottom = next ? next.getBoundingClientRect().top + window.scrollY - 10 : document.body.scrollHeight;
      const body = document.body.getBoundingClientRect();
      return { x: Math.max(body.left - 8, 0), y: top, width: Math.min(body.width + 16, 1000), height: Math.min(bottom - top, limit) };
    }, [heading, maxHeight, until]);
    if (!box) { console.log('missing', heading); continue; }
    await page.screenshot({ path: path.join(out, name + '.png'), clip: box, fullPage: true });
    console.log('saved', name, Math.round(box.height));
  }
  await browser.close();
})();
