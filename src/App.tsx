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
    'pletzer-lab-genome-browser:bacterium',
    bacteria[0],
  );

  const [viewState, setViewState] = useState<ViewModel>();
  const [conditionA, setConditionA] = useState<[number, number]>([0, 0]);
  const [conditionB, setConditionB] = useState<[number, number]>([1, 0]);

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

  // if selected bacterium does not exist, use default
  useEffect(() => {
    if (!myConf[bacterium]) setBacterium(bacteria[0]);
  }, [bacterium, bacteria, setBacterium]);

  // use default conditions on bacterium change
  useEffect(() => {
    setConditionA([0, 0]);
    setConditionB([1, 0]);
  }, [bacterium]);

  // full-refresh jbrowse required for bacterium or conditions changes
  useEffect(() => {
    const config = myConf[bacterium];
    if (!config) return;

    const state = myCreateViewState(
      buildConfig(config, {
        conditionA,
        conditionB,
        loc: [0, 5000], // TODO: keep loc constant for condition changes only
      }),
    );

    setViewState(state);
  }, [bacterium, conditionA, conditionB]);

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

  if (!viewState) return null;

  const coverage = myConf[bacterium].data.coverage;

  return (
    <>
      <header className="header">
        <div className="header-title">
          <img src="./pletzerlab-icon.webp" alt="Pletzer Lab Icon" />
          <h1>Pletzer Lab Genome Browser</h1>
        </div>

        <nav className="header-genome-chooser">
          <GenomeSelector
            bacteria={bacteria}
            selected={bacterium}
            onChange={setBacterium}
          />

          {coverage.length >= 2 && (
            <div className="header-condition-chooser">
              <ConditionSelect
                id="select-condition-a"
                label="Condition"
                coverage={coverage}
                value={conditionA}
                onChange={setConditionA}
              />
              <span>vs.</span>
              <ConditionSelect
                id="select-condition-b"
                label="Condition"
                coverage={coverage}
                value={conditionB}
                onChange={setConditionB}
              />
            </div>
          )}
        </nav>

        <div className="header-buttons-container">
          <div className="header-buttons">
            <BooleanCheckbox
              label="Log2"
              checked={logScaling}
              onChange={setLogScaling}
              hoverDescription="Use base-2 log scale"
            />
            <BooleanCheckbox
              label="Global"
              checked={globalScaling}
              onChange={setGlobalScaling}
              hoverDescription="Fix Y-axis between global minimum and maximum"
            />
            <BooleanCheckbox
              label="CDS+Amino"
              checked={colorByCds}
              onChange={setColorByCDS}
              hoverDescription="Colour by CDS and amino acids"
            />
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

function GenomeSelector({
  bacteria,
  selected,
  onChange,
}: {
  bacteria: string[];
  selected: string;
  onChange: (b: string) => void;
}) {
  return (
    <select value={selected} onChange={e => onChange(e.target.value)}>
      {bacteria.map(b => (
        <option key={b} value={b}>
          {b}
        </option>
      ))}
    </select>
  );
}

function ConditionSelect({
  label,
  coverage,
  value,
  onChange,
  id,
}: {
  label: string;
  coverage: string[][];
  value: [number, number];
  onChange: (val: [number, number]) => void;
  className?: string;
  id: string;
}) {
  return (
    <select
      className="condition-select"
      id={id}
      value={value.join(',')}
      onChange={e => {
        const [c, r] = e.target.value.split(',').map(Number);
        onChange([c, r]);
      }}
    >
      {coverage.map((arr, i) => (
        <optgroup key={`${id}--group-${i}`} label={`${label} ${i + 1}`}>
          {arr.map((fname, j) => (
            <option key={`${id}--group-${i}-item-${j}`} value={`${i},${j}`}>
              {fname}
            </option>
          ))}
        </optgroup>
      ))}
    </select>
  );
}

function BooleanCheckbox({
  label,
  checked,
  onChange,
  hoverDescription,
}: {
  label: string;
  checked: boolean;
  onChange: (b: boolean) => void;
  hoverDescription: string;
}) {
  const id = label.replace(/\s+/g, '');
  return (
    <div className="checkbox">
      <label htmlFor={id}>{label}</label>
      <input
        type="checkbox"
        id={id}
        checked={checked}
        onChange={e => onChange(e.target.checked)}
        title={hoverDescription}
      />
    </div>
  );
}
