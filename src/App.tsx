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

  // update view with selected dataset
  useEffect(() => {
    const config = myConf[bacterium];
    if (!config) {
      return;
    }
    const state = myCreateViewState(buildConfig(config));
    setViewState(state);
    resetCheckboxes();
  }, [bacterium]);

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
          <select onChange={event => setBacterium(event.target.value)}>
            {bacteria.map(b => (
              <option value={b} key={`dataset-${b}`}>
                {b}
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
                onChange={event => logScaleCheckboxHandler(event)}
              />
            </div>
            <div className="checkbox">
              <label htmlFor="aminoAcidsCheckbox">CDS+Amino</label>
              <input
                type="checkbox"
                id="aminoAcidsCheckbox"
                onChange={event => aminoAcidCheckboxHandler(event)}
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
    const checkboxes = document.querySelectorAll(
      '.header-buttons input[type="checkbox"]',
    ) as NodeListOf<HTMLInputElement>;

    checkboxes.forEach(box => {
      if (box.checked) {
        box.checked = false;
      }
    });
  }

  function aminoAcidCheckboxHandler(evt: React.ChangeEvent<HTMLInputElement>) {
    const checked = evt.target.checked;
    if (!viewState) return;

    const linearView = viewState.session.views.find(
      v => v.type === 'LinearGenomeView',
    );
    if (!linearView) return;
    linearView.setColorByCDS(checked);
  }

  function logScaleCheckboxHandler(evt: React.ChangeEvent<HTMLInputElement>) {
    const checked = evt.target.checked;
    if (!viewState) return;

    const linearView = viewState.session.views.find(
      v => v.type === 'LinearGenomeView',
    );
    if (!linearView) return;

    linearView.tracks.forEach(track => {
      /** @ts-expect-error display is 'any' -_- */
      track.displays.forEach(display => {
        if (display.type.includes('Wiggle')) {
          display.setScaleType(checked ? 'log' : 'linear');
        }
      });
    });
  }
}

export default App;
