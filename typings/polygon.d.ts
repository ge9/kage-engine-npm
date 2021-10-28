export declare class Polygon {
    get array(): ReadonlyArray<Readonly<{
        x: number;
        y: number;
        off: boolean;
    }>>;
    get length(): number;
    protected _array: {
        x: number;
        y: number;
        off: boolean;
    }[];
    constructor(param?: number | {
        x: number;
        y: number;
        off?: boolean;
    }[]);
    push(x: number, y: number, off?: boolean): void;
    set(index: number, x: number, y: number, off?: boolean): void;
    get(index: number): Readonly<{
        x: number;
        y: number;
        off: boolean;
    }>;
    reverse(): void;
    concat(poly: Polygon): void;
    shift(): void;
    unshift(x: number, y: number, off?: boolean): void;
    clone(): Polygon;
    translate(dx: number, dy: number): this;
    transformMatrix(a: number, b: number, c: number, d: number): this;
    /**
     * Scales by hypot(x, y) and rotates by atan2(y, x). Corresponds to multiplying x+yi on complex plane.
     */
    transformMatrix2(x: number, y: number): this;
    scale(factor: number): this;
    reflectX(): this;
    reflectY(): this;
    rotate90(): this;
    rotate180(): this;
    rotate270(): this;
    floor(): this;
}
