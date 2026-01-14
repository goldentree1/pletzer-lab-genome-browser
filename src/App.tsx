import type { ViewModel } from './types';
import { useState, useEffect } from 'react';
import { JBrowseLinearGenomeView } from '@jbrowse/react-linear-genome-view2';
// @ts-expect-error no types for font
import '@fontsource/roboto';
import myConf from './config';
import { buildConfig, myCreateViewState } from './jbrowse-custom';

function App() {
  const bacteria = Object.keys(myConf).sort();
  const [bacterium, setBacterium] = useState<string>(bacteria[0]);
  const [viewState, setViewState] = useState<ViewModel>();
  const [conditionA, setConditionA] = useState<[number, number]>([0, 0]);
  const [conditionB, setConditionB] = useState<[number, number]>([1, 0]);

  useEffect(() => {
    const config = myConf[bacterium];
    if (!config) return;

    const state = myCreateViewState(
      buildConfig(config, { conditionA, conditionB, loc: [0, 5000] }),
    );

    setViewState(state);
    resetCheckboxes();
  }, [bacterium, conditionA, conditionB]);

  if (!viewState) return null;

  return (
    <>
      <header className="header">
        <div className="header-title">
          <img src="./pletzerlab-icon.webp" alt="Pletzer Lab Icon" />
          <h1>Pletzer Lab Genome Browser</h1>
        </div>
        <nav className="header-genome-chooser">
          <select onChange={e => setBacterium(e.target.value)}>
            {bacteria.map(b => (
              <option key={b} value={b}>
                {b}
              </option>
            ))}
          </select>

          {myConf[bacterium].data.coverage.length >= 2 && (
            <div className="header-condition-chooser">
              <select
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

              <select
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
                onChange={logScaleCheckboxHandler}
              />
            </div>
            <div className="checkbox">
              <label htmlFor="aminoAcidsCheckbox">CDS+Amino</label>
              <input
                type="checkbox"
                id="aminoAcidsCheckbox"
                onChange={aminoAcidCheckboxHandler}
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

  function resetCheckboxes() {
    document
      .querySelectorAll<HTMLInputElement>(
        '.header-buttons input[type="checkbox"]',
      )
      .forEach(box => (box.checked = false));
  }

  function aminoAcidCheckboxHandler(evt: React.ChangeEvent<HTMLInputElement>) {
    const linearView = viewState?.session.views.find(
      v => v.type === 'LinearGenomeView',
    );
    if (!linearView) return;
    linearView.setColorByCDS(evt.target.checked);
  }

  function logScaleCheckboxHandler(evt: React.ChangeEvent<HTMLInputElement>) {
    const linearView = viewState?.session.views.find(
      v => v.type === 'LinearGenomeView',
    );
    if (!linearView) return;

    linearView.tracks.forEach(track => {
      /** @ts-expect-error display is 'any' */
      track.displays.forEach(display => {
        if (display.type.includes('Wiggle')) {
          display.setScaleType(evt.target.checked ? 'log' : 'linear');
        }
      });
    });
  }
}

export default App;
