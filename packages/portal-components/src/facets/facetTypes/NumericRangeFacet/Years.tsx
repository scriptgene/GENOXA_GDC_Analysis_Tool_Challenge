import React from "react";
import { getLowerAgeYears } from "../../utils";
import RangeInputWithPrefixedRanges from "./RangeInputWithPrefixedRanges";
import { NumericFacetProps } from "./types";

const Years: React.FC<NumericFacetProps> = ({
  field,
  valueLabel,
  rangeDatatype,
  hooks,
  clearValues,
  minimum = undefined,
  maximum = undefined,
  isFacetView,
  queryOptions,
  Chart,
}: NumericFacetProps) => {
  // minimum and maximum values if not undefined are in days so need to convert it to years
  const adjMinimum = minimum != undefined ? getLowerAgeYears(minimum) : 0;
  const adjMaximum = maximum != undefined ? getLowerAgeYears(maximum) : 89;
  const numBuckets = Math.round((adjMaximum - adjMinimum) / 10);

  return (
    <div className="flex flex-col w-100 space-y-2 px-1  mt-1 ">
      <RangeInputWithPrefixedRanges
        valueLabel={valueLabel}
        hooks={{ ...hooks }}
        rangeDatatype={rangeDatatype}
        units="years"
        minimum={adjMinimum}
        maximum={adjMaximum}
        numBuckets={numBuckets}
        field={field}
        clearValues={clearValues}
        isFacetView={isFacetView}
        queryOptions={queryOptions}
        Chart={Chart}
      />
    </div>
  );
};

export default Years;
