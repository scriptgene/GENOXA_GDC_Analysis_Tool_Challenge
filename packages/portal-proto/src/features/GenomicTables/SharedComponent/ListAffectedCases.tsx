import React, { useCallback, useMemo } from "react";
import { animated, useSpring } from "@react-spring/web";
import { useMeasure } from "react-use";
import { TableSubrowData } from "@gff/core";
import RatioWithSpring from "@/components/RatioWithSpring";
import { CircleIcon } from "@/utils/icons";

interface ListAffectedCasesProps {
  subData: ReadonlyArray<TableSubrowData>;
  subrowTitle: string;
}

interface TableSubrowDataWithRatio extends TableSubrowData {
  readonly ratio?: number;
}

const itemRatio = (item: TableSubrowData): number | undefined => {
  const { numerator, denominator } = item;
  const [n, d] = [
    denominator === 0 ? 0 : numerator,
    denominator === 0 ? 1 : denominator,
  ];

  if (denominator === 0) return undefined;
  return n / d;
};

const ListAffectedCases: React.FC<ListAffectedCasesProps> = ({
  subData,
  subrowTitle,
}: ListAffectedCasesProps) => {
  const [subRef, { width, height }] = useMeasure();

  const fudgeFactor = width / 60;

  const verticalSpring = useSpring({
    from: { opacity: 0, height: 0 },
    to: {
      opacity: 1,
      height: height + fudgeFactor,
    },
    immediate: true,
  });

  const subDataSorted = useMemo(
    () =>
      subData
        ?.map((x): TableSubrowDataWithRatio => ({ ...x, ratio: itemRatio(x) }))
        .sort((a: TableSubrowDataWithRatio, b: TableSubrowDataWithRatio) => {
          if (itemRatio(a) > itemRatio(b)) return -1;
          if (itemRatio(a) < itemRatio(b)) return 1;
          return 0;
        }),
    [subData],
  );

  const renderItems = useCallback(
    (item: TableSubrowDataWithRatio, index: number) => {
      return (
        <div className="flex items-center" key={`item-${index}`}>
          <CircleIcon size="0.65em" className="text-primary shrink-0 mr-2" />
          <RatioWithSpring
            item={{
              numerator: item.numerator ?? 0,
              denominator: item.denominator ?? 0,
              project: item.project,
            }}
            index={index}
          />
        </div>
      );
    },
    [],
  );

  return (
    <>
      <animated.div ref={subRef} className="absolute mt-2 ml-2 w-full">
        <div className="font-semibold text-[1rem] mb-2">{subrowTitle}</div>
        <div className="columns-4 font-content text-sm">
          {subDataSorted.map((item, i) => renderItems({ ...item }, i))}
        </div>
      </animated.div>
      <animated.div style={verticalSpring}></animated.div>
    </>
  );
};

export default ListAffectedCases;
