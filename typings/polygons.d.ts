import { Polygon } from "./polygon";
export declare class Polygons {
    array: Polygon[];
    constructor();
    clear(): void;
    push(polygon: Polygon): void;
    generateSVG(curve: boolean): string;
    generateEPS(): string;
}
