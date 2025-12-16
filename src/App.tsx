import { useState, useEffect } from 'react';
// @ts-expect-error no types for font
import '@fontsource/roboto';
import {
  createViewState,
  JBrowseLinearGenomeView,
} from '@jbrowse/react-linear-genome-view2';
type JBrowseConfig = Parameters<typeof createViewState>[0];

import config_NC011770 from './config-NC011770';
import config_NC008463 from './config-NC008463';

type ViewModel = ReturnType<typeof createViewState>;

function conf(datasetName: string): JBrowseConfig | null {
  if (datasetName === 'NC011770') {
    return config_NC011770;
  } else if (datasetName === 'NC008463') {
    return config_NC008463;
  } else {
    return null;
  }
}

function View() {
  const datasetNames = ['NC011770', 'NC008463'];
  const [configName, setConfigName] = useState(datasetNames[0]);
  const [viewState, setViewState] = useState<ViewModel>();

  useEffect(() => {
    const config = conf(configName);
    if (!config) {
      return;
    }
    const state = createViewState(config);
    setViewState(state);
  }, [configName]);

  if (!viewState) {
    return null;
  }

  return (
    <>
      <header className="header">
        <h1>Pletzer-JBrowse</h1>
        <div className="header-controls">
          <select onChange={evt => setConfigName(evt.target.value)}>
            {datasetNames.map(d => (
              <option value={d} key={`dataset-${d}`}>
                {d}
              </option>
            ))}
          </select>
          <div>
            <label htmlFor="log10Checkbox">Log10 Counts</label>
            <input type="checkbox" id="log10Checkbox" />
          </div>
          <div>
            <label htmlFor="sidebarCheckbox">Sidebar</label>
            <input type="checkbox" id="sidebarCheckbox" />
          </div>
          <button>Fullscreen</button>
        </div>
      </header>
      <div>
        <JBrowseLinearGenomeView viewState={viewState} />
      </div>
    </>
  );
}

export default View;
