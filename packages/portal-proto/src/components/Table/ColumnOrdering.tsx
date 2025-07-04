import { Dispatch, SetStateAction, useState } from "react";
import {
  Switch,
  ActionIcon,
  Tooltip,
  TextInput,
  Divider,
  Menu,
} from "@mantine/core";
import { Column, ColumnOrderState, Table } from "@tanstack/react-table";
import { isEqual } from "lodash";
import { humanify } from "@/utils/index";
import { NO_COLUMN_ORDERING_IDS } from "./utils";
import {
  DndContext,
  closestCenter,
  KeyboardSensor,
  PointerSensor,
  useSensor,
  useSensors,
  DragEndEvent,
} from "@dnd-kit/core";
import {
  arrayMove,
  SortableContext,
  sortableKeyboardCoordinates,
  verticalListSortingStrategy,
} from "@dnd-kit/sortable";
import { useSortable } from "@dnd-kit/sortable";
import { CSS } from "@dnd-kit/utilities";
import {
  restrictToParentElement,
  restrictToVerticalAxis,
  restrictToWindowEdges,
} from "@dnd-kit/modifiers";
import { DragIcon, ListIcon, SearchIcon, UndoIcon, XIcon } from "@/utils/icons";

function ColumnOrdering<TData>({
  table,
  handleColumnOrderingReset,
  columnOrder,
  setColumnOrder,
}: {
  table: Table<TData>;
  handleColumnOrderingReset: () => void;
  columnOrder: ColumnOrderState;
  setColumnOrder: Dispatch<SetStateAction<ColumnOrderState>>;
}): JSX.Element {
  const [showColumnMenu, setShowColumnMenu] = useState(false);
  const [searchValue, setSearchValue] = useState("");

  const isBackToDefaults =
    isEqual(table.initialState.columnOrder, columnOrder) &&
    isEqual(
      table.initialState.columnVisibility,
      Object.fromEntries(
        Object.entries(table.getState().columnVisibility).filter(
          ([, value]) => !value,
        ),
      ),
    );

  const sensors = useSensors(
    useSensor(PointerSensor, {
      activationConstraint: { distance: 5 },
    }),
    useSensor(KeyboardSensor, {
      coordinateGetter: sortableKeyboardCoordinates,
    }),
  );

  const handleDragEnd = (event: DragEndEvent) => {
    const { active, over } = event;

    if (active?.id && over?.id && active?.id !== over?.id) {
      setColumnOrder((items) => {
        const oldIndex = items.indexOf(active?.id as string);
        const newIndex = items.indexOf(over?.id as string);

        return arrayMove(items, oldIndex, newIndex);
      });
    }
  };

  return (
    <Menu
      position="bottom-end"
      middlewares={{ shift: false, flip: false }}
      onChange={setShowColumnMenu}
      offset={0}
      zIndex={400}
      withinPortal={true}
    >
      <Tooltip label="Customize Columns" disabled={showColumnMenu}>
        <Menu.Target>
          <ActionIcon
            variant="outline"
            size="lg"
            aria-label="Customize Columns"
            color="primary"
            data-testid="button-column-selector-box"
            className={`${
              showColumnMenu && "border-2 border-primary"
            } hover:bg-primary hover:text-base-max`}
          >
            {!showColumnMenu ? (
              <ListIcon size="1.5rem" aria-hidden="true" />
            ) : (
              <XIcon size="2rem" aria-hidden="true" />
            )}
          </ActionIcon>
        </Menu.Target>
      </Tooltip>

      <Menu.Dropdown
        className="bg-base-max border-2 border-primary p-2 rounded-md"
        data-testid="column-selector-popover-modal"
      >
        <div className="flex justify-between items-center">
          <span className="font-bold text-primary-darkest tracking-normal">
            Customize Columns
          </span>

          <Tooltip
            label="Restore defaults"
            disabled={isBackToDefaults}
            zIndex={400}
          >
            <span>
              <ActionIcon
                onClick={handleColumnOrderingReset}
                className={isBackToDefaults ? "invisible" : "border-0"}
                data-testid="restore-default-icon"
                aria-label="restore default column ordering"
              >
                <UndoIcon className="text-primary" size="1rem" />
              </ActionIcon>
            </span>
          </Tooltip>
        </div>
        <Divider color="base.2" className="mt-1" />
        <TextInput
          value={searchValue}
          onChange={(event) => setSearchValue(event.currentTarget.value.trim())}
          placeholder="Filter Columns"
          aria-label="Search input for columns"
          leftSection={<SearchIcon aria-hidden="true" />}
          className="mb-2 mt-4"
          data-testid="textbox-column-selector"
        />
        <DndContext
          sensors={sensors}
          collisionDetection={closestCenter}
          onDragEnd={handleDragEnd}
          modifiers={[
            restrictToVerticalAxis,
            restrictToWindowEdges,
            restrictToParentElement,
          ]}
        >
          <SortableContext
            items={table.getAllLeafColumns()}
            strategy={verticalListSortingStrategy}
          >
            <List
              columns={table.getAllLeafColumns()}
              searchValue={searchValue}
            />
          </SortableContext>
        </DndContext>
      </Menu.Dropdown>
    </Menu>
  );
}

function List<TData>({
  columns,
  searchValue,
}: {
  columns: Column<TData, unknown>[];
  searchValue: string;
}) {
  return (
    <ul>
      {columns
        .filter((column) => {
          if (!NO_COLUMN_ORDERING_IDS.includes(column.id)) {
            return humanify({ term: column.id })
              .toLowerCase()
              .includes(searchValue.toLowerCase());
          } else {
            return column.id;
          }
        })
        .map((column, index) => {
          return !NO_COLUMN_ORDERING_IDS.includes(column.id) ? (
            <DraggableColumnItem
              key={column.id}
              column={column}
              isNotLast={index < columns.length - 1}
            />
          ) : (
            <li className="hide" key={column.id} />
          );
        })}
    </ul>
  );
}

function DraggableColumnItem<TData>({
  column,
  isNotLast,
}: {
  column: Column<TData, unknown>;
  isNotLast: boolean;
}) {
  const { attributes, listeners, setNodeRef, transform, transition } =
    useSortable({ id: column.id });

  const style = {
    transform: CSS.Transform.toString(transform),
    transition,
  };

  return (
    <li
      ref={setNodeRef}
      style={style}
      {...listeners}
      className={` ${isNotLast ? "mb-2" : ""}`}
      data-testid={`column-selector-row-${column.id}`}
    >
      <div
        {...attributes}
        role={undefined}
        className="flex gap-2 items-center bg-nci-violet-lightest px-1 py-1.5 h-6 cursor-move"
      >
        <DragIcon size="1rem" className="text-primary" />
        <Switch
          labelPosition="left"
          label={
            typeof column?.columnDef?.header === "string"
              ? column.columnDef.header
              : humanify({ term: column.id })
          }
          classNames={{
            root: "w-full",
            body: "grow flex justify-between",
            label:
              "text-xs text-secondary-contrast-lighter font-medium tracking-normal",
          }}
          color="accent"
          checked={column.getIsVisible()}
          onChange={() => column.toggleVisibility()}
          onKeyDown={(e) => {
            if (e.code === "Space") {
              e.preventDefault();
              column.toggleVisibility();
            }
          }}
          size="xs"
          data-testid="switch-toggle"
        />
      </div>
    </li>
  );
}

export default ColumnOrdering;
