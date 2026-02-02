interface SelectProps {
  value: string;
  values: string[];
  onChange: (value: string) => void;
  id: string;
  className?: string;
  optionLabel?: (value: string) => string;
}

function Select({ value, values, onChange, id, className }: SelectProps) {
  return (
    <select
      className={className}
      id={id}
      value={value}
      onChange={e => onChange(e.target.value)}
    >
      {values.map((value, index) => (
        <option key={index} value={value}>
          {value}
        </option>
      ))}
    </select>
  );
}

export default Select;
