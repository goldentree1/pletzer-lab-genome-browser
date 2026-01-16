import type { ViewModel } from './types';
import { useState, useEffect } from 'react';
import { JBrowseLinearGenomeView } from '@jbrowse/react-linear-genome-view2';
// @ts-expect-error no types for font
import '@fontsource/roboto';
import myConf from './config';
import { buildConfig, myCreateViewState } from './jbrowse-custom';
import { useStoredStateBoolean, useStoredStateString } from './utils';

function App() {
  const bacteria: string[] = Object.keys(myConf).sort();
  const [bacterium, setBacterium] = useStoredStateString(
    'pletzer-genome-browser:bacterium',
    bacteria[0],
  );
  const [viewState, setViewState] = useState<ViewModel>();
  const [conditionA, setConditionA] = useState<[number, number]>([0, 0]);
  const [conditionB, setConditionB] = useState<[number, number]>([1, 0]);
  const [logScaling, setLogScaling] = useStoredStateBoolean(
    'pletzer-genome-browser:logScaling',
    true,
  );
  const [globalScaling, setGlobalScaling] = useStoredStateBoolean(
    'pletzer-genome-browser:globalScaling',
    true,
  );
  const [colorByCds, setColorByCDS] = useStoredStateBoolean(
    'pletzer-genome-browser:colorByCDS',
    false,
  );

  // config full reload (for things jbrowse cant handle via DOM)
  useEffect(() => {
    const config = myConf[bacterium];
    if (!config) return;

    const state = myCreateViewState(
      buildConfig(config, {
        conditionA,
        conditionB,
        loc: [0, 5000],
      }),
    );

    setViewState(state);
  }, [bacterium, conditionA, conditionB]);

  // edit browser via DOM manipulation (JBrowse functions)
  useEffect(() => {
    const linearView = viewState?.session.views.find(
      v => v.type === 'LinearGenomeView',
    );
    if (!linearView) return;

    // log scale + autoscaling
    linearView.tracks.forEach(track => {
      /** @ts-expect-error display is 'any' */
      track.displays.forEach(display => {
        if (display.type.includes('Wiggle')) {
          const wantedScale = logScaling ? 'log' : 'linear';
          if (display.scaleType !== wantedScale) {
            display.setScaleType(wantedScale);
          }
          if (display.autoscale !== (globalScaling ? 'global' : 'local')) {
            display.setAutoscale(globalScaling ? 'global' : 'local');
          }
        }
      });
    });

    // colour by CDS
    linearView.setColorByCDS(colorByCds);
  }, [logScaling, colorByCds, globalScaling, viewState]);

  useEffect(() => {
    // revert conditions to defaults on bacterium change
    setConditionA([0, 0]);
    setConditionB([1, 0]);
  }, [bacterium]);

  if (!viewState) return null;

  return (
    <>
      <header className="header">
        <div className="header-title">
          <img src="./pletzerlab-icon.webp" alt="Pletzer Lab Icon" />
          <h1>Pletzer Lab Genome Browser</h1>
        </div>
        <nav className="header-genome-chooser">
          <select
            defaultValue={bacteria.find(b => b === bacterium) || bacteria[0]}
            onChange={e => setBacterium(e.target.value)}
          >
            {bacteria.map(b => (
              <option key={b} value={b}>
                {b}
              </option>
            ))}
          </select>

          {myConf[bacterium].data.coverage.length >= 2 && (
            <div className="header-condition-chooser">
              <select
                className="condition-select"
                id="select-condition-a"
                onChange={e => {
                  const [c, r] = e.target.value.split(',').map(Number);
                  setConditionA([c, r]);
                }}
                value={conditionA.join(',')}
              >
                {myConf[bacterium].data.coverage.map((arr, i) => (
                  <optgroup
                    key={`condA-group-${i}`}
                    label={`Condition ${i + 1}`}
                  >
                    {arr.map((fname, j) => (
                      <option key={`condA-${i}-${j}`} value={`${i},${j}`}>
                        {fname}
                      </option>
                    ))}
                  </optgroup>
                ))}
              </select>
              <span>vs.</span>
              <select
                className="condition-select"
                id="select-condition-b"
                onChange={e => {
                  const [c, r] = e.target.value.split(',').map(Number);
                  setConditionB([c, r]);
                }}
                value={conditionB.join(',')}
              >
                {myConf[bacterium].data.coverage.map((arr, i) => (
                  <optgroup
                    key={`condB-group-${i}`}
                    label={`Condition ${i + 1}`}
                  >
                    {arr.map((fname, j) => (
                      <option key={`condB-${i}-${j}`} value={`${i},${j}`}>
                        {fname}
                      </option>
                    ))}
                  </optgroup>
                ))}
              </select>
            </div>
          )}
        </nav>

        <div className="header-buttons-container">
          <div className="header-buttons">
            <div className="checkbox">
              <label htmlFor="logScaleCheckbox">Log Scale</label>
              <input
                type="checkbox"
                id="logScaleCheckbox"
                checked={logScaling}
                onChange={e => setLogScaling(e.target.checked)}
              />
            </div>
            <div className="checkbox">
              <label htmlFor="globalScaleCheckbox">Global Scaling</label>
              <input
                type="checkbox"
                id="globalScaleCheckbox"
                checked={globalScaling}
                onChange={e => setGlobalScaling(e.target.checked)}
              />
            </div>

            <div className="checkbox">
              <label htmlFor="aminoAcidsCheckbox">CDS+Amino</label>
              <input
                type="checkbox"
                id="aminoAcidsCheckbox"
                checked={colorByCds}
                onChange={e => setColorByCDS(e.target.checked)}
              />
            </div>
          </div>
        </div>
      </header>

      <div className="jbrowse-container">
        {viewState && <JBrowseLinearGenomeView viewState={viewState} />}
      </div>
    </>
  );
}

export default App;
