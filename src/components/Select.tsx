interface SelectProps {
  label: string;
  value: string;
  setValue: (value: string) => void;
  options: string[];
}

function Select({ label, value, setValue, options }: SelectProps) {
  return (
    <div className="select">
      <label htmlFor="select">{label}</label>
      <select value={value} onChange={e => setValue(e.target.value)}>
        {options.map((option, index) => (
          <option key={index} value={option}>
            {option}
          </option>
        ))}
      </select>
    </div>
  );
}

export default Select;
