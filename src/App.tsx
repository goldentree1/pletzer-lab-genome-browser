import { useState, useEffect } from 'react';
// @ts-expect-error no types for font
import '@fontsource/roboto';
import {
  createViewState,
  JBrowseLinearGenomeView,
} from '@jbrowse/react-linear-genome-view2';

import { config } from './config_NC011770';

type ViewModel = ReturnType<typeof createViewState>;

function View() {
  const [viewState, setViewState] = useState<ViewModel>();

  useEffect(() => {
    const state = createViewState(config);

    setViewState(state);
  }, []);

  if (!viewState) {
    return null;
  }

  return (
    <>
      <h1>JBrowse 2 React Linear Genome View Demo w/ vite</h1>
      <JBrowseLinearGenomeView viewState={viewState} />
    </>
  );
}

export default View;
