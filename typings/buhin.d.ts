export declare class Buhin {
    protected hash: {
        [name: string]: string;
    };
    constructor();
    set(name: string, data: string): void;
    search(name: string): string;
    push(name: string, data: string): void;
}
