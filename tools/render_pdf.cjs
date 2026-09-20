/* Optional PDF renderer. Requires Playwright and a running generated-docs server. */
const fs = require('fs');
const path = require('path');

async function render({browser, baseUrl, output}) {
  const page = await browser.newPage({viewport:{width:1200,height:1000}});
  const errors=[];
  page.on('pageerror', e=>errors.push(e.message));
  page.on('response', r=>{if(r.status()>=400) errors.push(r.status()+' '+r.url());});
  await page.goto(baseUrl+'/paper/manuscript.html',{waitUntil:'networkidle'});
  await page.waitForFunction(()=>window.MathJax && MathJax.startup && MathJax.startup.promise);
  await page.evaluate(()=>MathJax.startup.promise);
  const counts=await page.evaluate(()=>({
    math:document.querySelectorAll('mjx-container').length,
    mathErrors:document.querySelectorAll('mjx-merror,[data-mml-node="merror"]').length,
    tables:document.querySelectorAll('table').length,
    placeholders:document.body.innerText.includes('MATHPLACEHOLDER')
  }));
  if(errors.length || counts.mathErrors || counts.placeholders || counts.math<1900 || counts.tables!==28)
    throw new Error(JSON.stringify({errors,counts}));
  await page.emulateMedia({media:'print'});
  await page.evaluate(async()=>{await document.fonts.ready;await new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve)));});
  await page.pdf({path:output,format:'A4',printBackground:true,preferCSSPageSize:true,
    displayHeaderFooter:true,outline:true,tagged:true,
    headerTemplate:'<div style="font-size:8px;width:100%;text-align:center;color:#666">FrailtyBLUP · Manuscript 0.12</div>',
    footerTemplate:'<div style="font-size:9px;width:100%;text-align:center;color:#666"><span class="pageNumber"></span> / <span class="totalPages"></span></div>'});
  await page.emulateMedia({media:'screen'});
  const pages=[];
  for(const route of ['/index.html','/docs/data.html','/docs/api.html','/docs/examples.html','/paper/manuscript.html']) {
    await page.setViewportSize({width:390,height:844});
    await page.goto(baseUrl+route,{waitUntil:'networkidle'});
    await page.evaluate(()=>MathJax.startup.promise);
    const size=await page.evaluate(()=>({width:innerWidth,documentWidth:document.documentElement.scrollWidth}));
    if(size.documentWidth>size.width+1) throw new Error('Mobile overflow '+route+' '+JSON.stringify(size));
    pages.push({route,...size});
  }
  await page.close();
  const report={passed:errors.length===0,counts,pages,errors};
  fs.writeFileSync(path.join(path.dirname(output),'render_checks.json'),JSON.stringify(report,null,2)+'\n');
  console.log(JSON.stringify(report,null,2));
  return report;
}

if(require.main===module) {
  (async()=>{
    const {chromium}=require('playwright');
    const browser=await chromium.launch();
    try {await render({browser,baseUrl:process.argv[2] || 'http://127.0.0.1:8000',output:process.argv[3] || 'paper/paper.pdf'});}
    finally {await browser.close();}
  })().catch(e=>{console.error(e);process.exit(1);});
}
module.exports={render};
