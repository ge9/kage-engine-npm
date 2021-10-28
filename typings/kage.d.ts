import { Buhin } from "./buhin";
import { Polygons } from "./polygons";
import { stretch, Stroke } from "./stroke";
import { Font } from "./font";
export declare enum KShotai {
    kMincho = 0,
    kGothic = 1
}
export declare class Kage {
    static Buhin: typeof Buhin;
    static Polygons: typeof Polygons;
    kMincho: KShotai;
    kGothic: KShotai;
    kFont: Font;
    get kShotai(): KShotai;
    set kShotai(shotai: KShotai);
    get kUseCurve(): boolean;
    set kUseCurve(value: boolean);
    kBuhin: Buhin;
    stretch: typeof stretch;
    constructor(size?: number);
    makeGlyph(polygons: Polygons, buhin: string): void;
    makeGlyph2(polygons: Polygons, data: string): void;
    makeGlyph3(data: string): Polygons[];
    makeGlyphSeparated(data: string[]): Polygons[];
    protected getEachStrokes(glyphData: string): Stroke[];
    protected getEachStrokesOfBuhin(buhin: string, x1: number, y1: number, x2: number, y2: number, sx: number, sy: number, sx2: number, sy2: number): Stroke[];
    protected getBox(strokes: Stroke[]): {
        minX: number;
        maxX: number;
        minY: number;
        maxY: number;
    };
}
