import type { ViewModel } from './types';
import { useState, useEffect, useRef } from 'react';
import { JBrowseLinearGenomeView } from '@jbrowse/react-linear-genome-view2';
// @ts-expect-error no types for font
import '@fontsource/roboto';
import myConf from './config';
import { buildConfig, myCreateViewState } from './jbrowse-custom-hybrid';
import { useStoredStateBoolean, useStoredStateString } from './utils';
import Checkbox from './components/Checkbox';
import Select from './components/Select';
import { autorun } from 'mobx';
import ConditionsSelect from './components/ConditionSelect';

function App() {
  const genomes: string[] = Object.keys(myConf).sort();
  const [genome, setGenome] = useStoredStateString(
    'pletzer-lab-genome-browser:genome',
    genomes[0],
  );

  const [viewState, setViewState] = useState<ViewModel>();
  const [conditionA, setConditionA] = useState<[number, number]>([0, 0]);
  const [conditionB, setConditionB] = useState<[number, number]>([1, 0]);
  const [norms, setNorms] = useState(['none']);
  const loc = useRef<string | null>(null);

  const [preferredExperiment, setPreferredExperiment] = useStoredStateString(
    'pletzer-lab-genome-browser:preferredExperiment',
    '',
  );
  const [genesLabelType, setGenesLabelType] = useStoredStateString(
    'pletzer-lab-genome-browser:genesLabelType',
    'name',
  );
  const [normType, setNormType] = useStoredStateString(
    'pletzer-lab-genome-browser:normType',
    'none',
  );
  const [logScaling, setLogScaling] = useStoredStateBoolean(
    'pletzer-lab-genome-browser:logScaling',
    true,
  );
  const [globalScaling, setGlobalScaling] = useStoredStateBoolean(
    'pletzer-lab-genome-browser:globalScaling',
    true,
  );
  const [colorByCds, setColorByCDS] = useStoredStateBoolean(
    'pletzer-lab-genome-browser:colorByCDS',
    false,
  );

  const experiments = Object.keys(myConf[genome].data.experiments);

  const experiment = experiments.includes(preferredExperiment)
    ? preferredExperiment
    : experiments[0];
  const activeExperiment =
    experiment && experiments.includes(experiment)
      ? experiment
      : experiments[0];

  const coverage = activeExperiment
    ? myConf[genome].data.experiments[activeExperiment].coverage
    : [];

  const coverageConditionNames = activeExperiment
    ? myConf[genome].data.experiments[activeExperiment].coverage_condition_names
    : [];

  useEffect(() => {
    const linearView = viewState?.session.views.find(
      v => v.type === 'LinearGenomeView',
    );
    if (!linearView) return;

    const dispose = autorun(() => {
      loc.current = getLocString(linearView) || null;
    });

    return () => dispose();
  }, [viewState]);

  // DOM manipulation works for some JBrowse features - like scaling and CDS colouring.
  // Smoother and faster than a full-refresh, so use it where possible.
  useEffect(() => {
    const linearView = viewState?.session.views.find(
      v => v.type === 'LinearGenomeView',
    );
    if (!linearView) return;

    linearView.tracks.forEach(track => {
      /** @ts-expect-error display is 'any' */
      track.displays.forEach(display => {
        if (display.type.includes('Wiggle')) {
          const wantedScale = logScaling ? 'log' : 'linear';
          if (display.scaleType !== wantedScale)
            display.setScaleType(wantedScale);
          if (display.autoscale !== (globalScaling ? 'global' : 'local'))
            display.setAutoscale(globalScaling ? 'global' : 'local');
        }
      });
    });

    linearView.setColorByCDS(colorByCds);
  }, [viewState, logScaling, globalScaling, colorByCds]);

  // revert default conditions + loc on genome change
  useEffect(() => {
    setConditionA([0, 0]);
    setConditionB([1, 0]);
    loc.current = null;

    // const exps = Object.keys(myConf[genome].data.experiments);
    // setExperiment(exps[0] ?? '');
  }, [genome]);

  // full-refresh jbrowse required for genome or conditions changes
  useEffect(() => {
    const newLoc = [0, 5000];
    if (loc.current) {
      const [newStartLoc, newEndLoc] = loc.current
        .split(':')[1]
        .split('..')
        .map(Number);
      newLoc[0] = newStartLoc;
      newLoc[1] = newEndLoc;
      loc.current = null;
    }

    const config = myConf[genome];
    if (!config) {
      setGenome(genomes[0]);
      return;
    }

    setNorms(config.norms);

    if (!config.norms.includes(normType)) {
      setNormType(config.norms[0]);
      return;
    }

    // console.log('using loc:', newLoc);

    const state = myCreateViewState(
      buildConfig(config, {
        experiment,
        conditionA,
        conditionB,
        loc: newLoc,
        normType,
        genesLabelType,
      }),
    );

    setViewState(state);
  }, [genome, experiment, conditionA, conditionB, normType, genesLabelType]); // eslint-disable-line react-hooks/exhaustive-deps

  return (
    <>
      <header className="header">
        <nav className="header-genome-chooser">
          <div className="header-title">
            <img src="./pletzerlab-icon.webp" alt="Pletzer Lab Icon" />
            {/*<h1>Pletzer Lab Genome Browser</h1>*/}
          </div>
          <Select
            id="genome-select"
            values={genomes}
            value={genome}
            onChange={setGenome}
          />
          {experiments.length >= 1 && (
            <select
              className="experiment-select"
              value={experiment}
              onChange={e => setPreferredExperiment(e.target.value)}
            >
              {experiments.map(exp => (
                <option key={exp} value={exp}>
                  {exp}
                </option>
              ))}
            </select>
          )}

          {coverage.length >= 1 && (
            <div className="header-condition-chooser">
              <ConditionsSelect
                id="select-condition-a"
                label="Condition"
                coverage={coverage}
                coverageConditionNames={coverageConditionNames}
                value={conditionA}
                onChange={setConditionA}
              />
              {coverage.length > 1 && <span>vs.</span>}
              {coverage.length > 1 && (
                <ConditionsSelect
                  id="select-condition-b"
                  label="Condition"
                  coverage={coverage}
                  coverageConditionNames={coverageConditionNames}
                  value={conditionB}
                  onChange={setConditionB}
                />
              )}
            </div>
          )}
        </nav>

        <div className="header-buttons-container">
          <div className="header-buttons">
            <div>
              <label htmlFor="genes-label-type-select">Gene Labels:</label>
              <Select
                id="genes-label-type-select"
                className="genes-label-type-select"
                values={['name', 'locus_tag', 'old_locus_tag']}
                value={genesLabelType}
                onChange={setGenesLabelType}
              />
            </div>
            <div>
              <label htmlFor="norm-type-select">Normalisation: </label>
              <Select
                id="norm-type-select"
                className="norm-type-select"
                values={norms}
                value={normType}
                onChange={setNormType}
              />
            </div>
            <Checkbox
              label="Log2"
              checked={logScaling}
              onChange={setLogScaling}
              hoverDescription="Use base-2 log scale"
            />
            <Checkbox
              label="Global Y-axis"
              checked={globalScaling}
              onChange={setGlobalScaling}
              hoverDescription="Fix Y-axis between global minimum and maximum values"
            />
            <Checkbox
              label="CDS+Amino"
              checked={colorByCds}
              onChange={setColorByCDS}
              hoverDescription="Colour by CDS and amino acids (shows codon alignments by colour)"
            />
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

// eslint-disable-next-line @typescript-eslint/no-explicit-any
function getLocString(linearView: any) {
  if (!linearView.displayedRegions.length) return;
  const region = linearView.displayedRegions[0];

  const bpPerPx = linearView.bpPerPx;
  const offsetPx = linearView.offsetPx;
  const widthPx = linearView.width;

  const start = region.start + Math.floor(offsetPx * bpPerPx) + 1;
  const end = region.start + Math.floor((offsetPx + widthPx) * bpPerPx);

  const s = Math.max(start, 1);
  const e = Math.min(end, region.end);

  return region.reversed
    ? `${region.refName}:${e}..${s}`
    : `${region.refName}:${s}..${e}`;
}
