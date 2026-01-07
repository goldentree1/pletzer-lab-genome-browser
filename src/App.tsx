import { useState, useEffect } from 'react';
// @ts-expect-error no types for font
import '@fontsource/roboto';
import { JBrowseLinearGenomeView } from '@jbrowse/react-linear-genome-view2';
import type { JBrowseConfig, ViewModel } from './types';

import config_LESB58 from './config-LESB58';
import config_GCF_000014625 from './config-GCF_000014625';
import config_GCF_000006765_1 from './config-GCF_000006765.1';
import config_GCF_000013425_1 from './config-GCF_000013425.1';
import config_GCF_000013465_1 from './config-GCF_000013465.1';
import config_GCF_000281535_2 from './config-GCF_000281535.2';
import config_GCF_031932345_1 from './config-GCF_031932345.1';
import { myCreateViewState } from './makeConfig';

function App() {
  const datasetNames = [
    'P.aeruginosa LESB58',
    'P.aeruginosa PA14',
    'P.aeruginosa PAO1',
    'S.aureus HG001 (NCTC 8325)',
    'S.aureus USA300LAC',
    'K.pneumoniae KPNIH1',
    'A.baumannii AB5075-UW',
  ];
  const [configName, setConfigName] = useState(datasetNames[0]);
  const [viewState, setViewState] = useState<ViewModel>();

  // update view with selected dataset
  useEffect(() => {
    const config = getConf(configName);
    if (!config) {
      return;
    }
    const state = myCreateViewState(config);
    setViewState(state);

    /** TEMPORARY DISGUSTING
     * PIECE OF CRAP CODE - just undoes all checkboxes
     * on conf change (we should instead be APPLYING checkbox to the conf somehow)
     */
    const checboxesTEMPORARY = document.querySelectorAll(
      '.header-buttons input[type="checkbox"]',
    );
    checboxesTEMPORARY.forEach((el: Element) => {
      const el2 = el as HTMLInputElement;
      if (el2.checked) {
        el2.checked = false;
      }
    });
  }, [configName]);

  if (!viewState) {
    return null;
  }

  return (
    <>
      <header className="header">
        <div className="header-title">
          <img src="./pletzerlab-icon.webp" alt="Pletzer Lab Icon" />
          <h1>Pletzer Lab Genome Browser</h1>
        </div>
        <nav className="header-genome-chooser">
          <select onChange={event => setConfigName(event.target.value)}>
            {datasetNames.map(datasetName => (
              <option value={datasetName} key={`dataset-${datasetName}`}>
                {datasetName}
              </option>
            ))}
          </select>
        </nav>
        <div className="header-buttons-container">
          <div className="header-buttons">
            <div className="checkbox">
              <label htmlFor="logScaleCheckbox">Log Scale</label>
              <input
                type="checkbox"
                id="logScaleCheckbox"
                onChange={evt => {
                  const checked = evt.target.checked;
                  if (!viewState) return;

                  const linearView = viewState.session.views.find(
                    v => v.type === 'LinearGenomeView',
                  );
                  if (!linearView) return;

                  linearView.tracks.forEach(track => {
                    /** @ts-expect-error display is 'any' -_- */
                    track.displays.forEach(display => {
                      if (isWiggleDisplay(display)) {
                        display.setScaleType(checked ? 'log' : 'linear');
                      }
                    });
                  });
                }}
              />
            </div>
            <div className="checkbox">
              <label htmlFor="aminoAcidsCheckbox">CDS+Amino</label>
              <input
                type="checkbox"
                id="aminoAcidsCheckbox"
                onChange={evt => {
                  const checked = evt.target.checked;
                  if (!viewState) return;

                  const linearView = viewState.session.views.find(
                    v => v.type === 'LinearGenomeView',
                  );
                  if (!linearView) return;
                  linearView.setColorByCDS(checked);
                }}
              />
            </div>
          </div>
        </div>
      </header>
      <div className="jbrowse-container">
        <JBrowseLinearGenomeView viewState={viewState} />
      </div>
    </>
  );
}

export default App;

function getConf(datasetName: string): JBrowseConfig | null {
  if (datasetName === 'P.aeruginosa LESB58') {
    return config_LESB58;
  } else if (datasetName === 'P.aeruginosa PA14') {
    return config_GCF_000014625;
  } else if (datasetName === 'P.aeruginosa PAO1') {
    return config_GCF_000006765_1;
  } else if (datasetName === 'S.aureus HG001 (NCTC 8325)') {
    return config_GCF_000013425_1;
  } else if (datasetName === 'S.aureus USA300LAC') {
    return config_GCF_000013465_1;
  } else if (datasetName === 'K.pneumoniae KPNIH1') {
    return config_GCF_000281535_2;
  } else if (datasetName === 'A.baumannii AB5075-UW') {
    return config_GCF_031932345_1;
  } else {
    return null;
  }
}

function isWiggleDisplay(d: { type: string }) {
  return d.type.includes('Wiggle');
}
