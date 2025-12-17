import { useState, useEffect } from 'react';
// @ts-expect-error no types for font
import '@fontsource/roboto';
import {
  createViewState,
  JBrowseLinearGenomeView,
} from '@jbrowse/react-linear-genome-view2';
import type { JBrowseConfig, ViewModel } from './types';

import config_LESB58 from './config-LESB58';
import config_GCF_000014625 from './config-GCF_000014625';
import config_GCF_000006765_1 from './config-GCF_000006765.1';
import config_GCF_000013425_1 from './config-GCF_000013425.1';
import config_GCF_000013465_1 from './config-GCF_000013465.1';
import config_GCF_000281535_2 from './config-GCF_000281535.2';
import config_GCF_031932345_1 from './config-GCF_031932345.1';
import { makeConfig } from './makeConfig';

function App() {
  const datasetNames = [
    'LESB58',
    'GCF_000014625',
    'GCF_000006765.1',
    'GCF_000013425.1',
    'GCF_000013465.1',
    'GCF_000281535.2',
    'GCF_031932345.1',
  ];
  const [configName, setConfigName] = useState(datasetNames[0]);
  const [viewState, setViewState] = useState<ViewModel>();

  // update view with selected dataset
  useEffect(() => {
    const config = getConf(configName);
    if (!config) {
      return;
    }

    const state = createViewState(makeConfig(config));
    setViewState(state);
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
        <div className="header-genome-chooser">
          <select onChange={event => setConfigName(event.target.value)}>
            {datasetNames.map(datasetName => (
              <option value={datasetName} key={`dataset-${datasetName}`}>
                {datasetName}
              </option>
            ))}
          </select>
        </div>
        <div className="header-buttons-container">
          <div className="header-buttons">
            <div>
              <label htmlFor="logScaleCheckbox">Log Scale</label>
              <input type="checkbox" id="logScaleCheckbox" />
            </div>
            <div>
              <label htmlFor="legendCheckbox">Legend</label>
              <input type="checkbox" id="legendCheckbox" />
            </div>
            <div>
              <button>Fullscreen</button>
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
  if (datasetName === 'LESB58') {
    return config_LESB58;
  } else if (datasetName === 'GCF_000014625') {
    return config_GCF_000014625;
  } else if (datasetName === 'GCF_000006765.1') {
    return config_GCF_000006765_1;
  } else if (datasetName === 'GCF_000013425.1') {
    return config_GCF_000013425_1;
  } else if (datasetName === 'GCF_000013465.1') {
    return config_GCF_000013465_1;
  } else if (datasetName === 'GCF_000281535.2') {
    return config_GCF_000281535_2;
  } else if (datasetName === 'GCF_031932345.1') {
    return config_GCF_031932345_1;
  } else {
    return null;
  }
}
