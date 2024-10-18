// Description: 用来判断点在路径上 对应path的各种图形
module PointOnGraph {
    export class Point {
        x: number = 0;
        y: number = 0;
    }

    /**
     * 判断点是否在直线上
     * @param px 点的x坐标
     * @param py 点的y坐标
     * @param lx1 直线的起点x坐标
     * @param ly1 直线的起点y坐标
     * @param lx2 直线的终点x坐标
     * @param ly2 直线的终点y坐标
     * @param tolerance 误差范围
     * @returns {Point | undefined} 返回点在直线上的坐标(最近点)
     */
    export function pointOnLine(px: number, py: number, lx1: number, ly1: number, lx2: number, ly2: number, tolerance: number = 0): Point | undefined {
        // 点到直线的距离公式 |Ax + By + C| / sqrt(A^2 + B^2)
        // 其中
        // A 为 点1到点2的y差值
        // B 为 点1到点2的x差值
        // C 为 点1到点2的叉乘
        const A = ly2 - ly1;
        const B = lx1 - lx2;
        const C = lx2 * ly1 - lx1 * ly2;
        const dis = Math.abs(A * px + B * py + C) / Math.sqrt(A * A + B * B);
        if (dis <= tolerance) {
            // 点到直线的距离在误差范围内
            // 求交点 x = (B * B * px - A * B * py - A * C) / (A * A + B * B) y = (A * A * py - A * B * px - B * C) / (A * A + B * B)
            const x = (B * B * px - A * B * py - A * C) / (A * A + B * B);
            const y = (A * A * py - A * B * px - B * C) / (A * A + B * B);
            if (x >= Math.min(lx1, lx2) && x <= Math.max(lx1, lx2) && y >= Math.min(ly1, ly2) && y <= Math.max(ly1, ly2)) {
                return {x: x, y: y};
            }
        }
        return undefined;
    }

    /**
     * 判断点是否在矩形线条上
     * @param px 点的x坐标
     * @param py 点的y坐标
     * @param x 矩形的起点x坐标
     * @param y 矩形的起点y坐标
     * @param width 矩形的宽度
     * @param height 矩形的高度
     * @param tolerance 误差范围
     * @returns {Point | undefined} 点是否在矩形线条上
     */
    export function pointOnRect(px: number, py: number, x: number, y: number, width: number, height: number, tolerance: number = 0): Point | undefined {
        // 贴近上边 判断x在上边的范围内
        if (Math.abs(py - y) <= tolerance && px >= x && px <= x + width) {
            return {x: px, y: y};
        }
        // 贴近下边 判断x在下边的范围内
        if (Math.abs(py - (y + height)) <= tolerance && px >= x && px <= x + width) {
            return {x: px, y: y + height};
        }
        // 贴近左边 判断y在左边的范围内
        if (Math.abs(px - x) <= tolerance && py >= y && py <= y + height) {
            return {x: x, y: py};
        }
        // 贴近右边 判断y在右边的范围内
        if (Math.abs(px - (x + width)) <= tolerance && py >= y && py <= y + height) {
            return {x: x + width, y: py};
        }
        return undefined;
    }

    /**
     * 判断点是否在多条线段上
     * @param px 点的x坐标
     * @param py 点的y坐标
     * @param points 线段的点集
     * @param tolerance 误差范围
     * @returns {Point | undefined} 点是否在多条线段上
     */
    export function pointOnPolyline(px: number, py: number, points: number[], tolerance: number = 0): Point | undefined {
        for (let i = 0; i < points.length - 2; i += 2) {
            const point = pointOnLine(px, py, points[i], points[i + 1], points[i + 2], points[i + 3], tolerance);
            if (point) {
                return point;
            }
        }
        return undefined;
    }

    /**
     * 判断点是否在多边形线条上
     * @param px 点的x坐标
     * @param py 点的y坐标
     * @param points 多边形的点集
     * @param tolerance 误差范围
     * @returns {Point | undefined} 点是否在多边形线条上
     */
    export function pointOnPolygon(px: number, py: number, points: number[], tolerance: number = 0): Point | undefined {
        for (let i = 0; i < points.length - 2; i += 2) {
            const point = pointOnLine(px, py, points[i], points[i + 1], points[i + 2], points[i + 3], tolerance);
            if (point) {
                return point;
            }
        }
        const point = pointOnLine(px, py, points[points.length - 2], points[points.length - 1], points[0], points[1], tolerance);
        if (point) {
            return point;
        }
        return undefined;
    }

    /**
     * 判断点是否在扇形线条上
     * @param px 点的x坐标
     * @param py 点的y坐标
     * @param cx 圆的中心点x坐标
     * @param cy 圆的中心点y坐标
     * @param rx 扇形的x轴半径
     * @param ry 扇形的y轴半径
     * @param clockwise 扇形的方向
     * @param startAngle 扇形的起始角度
     * @param endAngle 扇形的结束角度
     * @param withStartLine 是否包含起点到中心点的线
     * @param withEndLine 是否包含终点到中心点的线
     * @param tolerance 误差范围
     * @returns {Point | undefined} 点是否在扇形线条上
     */
    export function pointOnEllipse(px: number, py: number, cx: number, cy: number, rx: number, ry: number, clockwise: boolean, startAngle: number, endAngle: number, withStartLine: boolean, withEndLine: boolean, tolerance: number = 0): Point | undefined {
        const startAngleRadians = startAngle * Math.PI / 180;
        const endAngleRadians = endAngle * Math.PI / 180;
        // 先判断点是否在那两条线上
        if (withStartLine) {
            const point = pointOnLine(px, py, cx, cy, cx + rx * Math.cos(startAngleRadians), cy + ry * Math.sin(startAngleRadians), tolerance);
            if (point) {
                return point;
            }
        }
        if (withEndLine) {
            const point = pointOnLine(px, py, cx, cy, cx + rx * Math.cos(endAngleRadians), cy + ry * Math.sin(endAngleRadians), tolerance);
            if (point) {
                return point;
            }
        }
        // 判断点是否在椭圆扇形上
        // 用过圆心与输入的点的连线判断点是否在椭圆上
        const angle = Math.atan2(py - cy, px - cx);
        // 判断角度是否在起始角度和结束角度之间
        if (clockwise) {
            if (angle < startAngleRadians || angle > endAngleRadians) {
                return undefined;
            }
        } else {
            if (angle > startAngleRadians || angle < endAngleRadians) {
                return undefined;
            }
        }
        const angleTan = Math.tan(angle);
        const angleTan2 = angleTan * angleTan;
        // 在范围内的话 先用角度计算椭圆上的点
        // x = - a * b / sqrt(b^2 + a^2 * tan^2(angle))
        // y = - a * b * tan(angle) / sqrt(b^2 + a^2 * tan^2(angle))
        const x = cx - rx * ry / Math.sqrt(ry * ry + rx * rx * angleTan2);
        const y = cy - rx * ry * angleTan / Math.sqrt(ry * ry + rx * rx * angleTan2);
        // 判断 输入点到 该点的距离是否在误差范围内
        const dis = Math.sqrt((px - x) * (px - x) + (py - y) * (py - y));
        if (dis <= tolerance) {
            return {x: x, y: y};
        }
        return undefined;
    }

    /**
     * 聚和 对应 C(n, m)
     * @param n
     * @param m
     * @private
     */
    function combination(n: number, m: number): number {
        // 组合数 C(n, m) = n! / m! * (n - m)!
        let result = 1;
        for (let i = 1; i <= m; i++) {
            result *= (n - i + 1) / i;
        }
        return result;
    }

    /**
     * 计算贝塞尔曲线上的点
     * @param t 贝塞尔曲线的参数
     * @param points 贝塞尔曲线的点集
     * @param n 贝塞尔曲线的阶数
     * @private
     */
    function bezierCurve(t: number, points: number[], n: number): Point {
        // n 阶贝塞尔曲线的公式 B(t) = Σ C(n, i) * (1 - t)^(n - i) * t^i * P(i)
        // x = Σ C(n, i) * (1 - t)^(n - i) * t^i * P(i).x
        // y = Σ C(n, i) * (1 - t)^(n - i) * t^i * P(i).y
        let x = 0;
        let y = 0;
        for (let i = 0; i <= n; i++) {
            const c = combination(n, i);
            const p = Math.pow(1 - t, n - i) * Math.pow(t, i);
            x += c * p * points[i * 2];
            y += c * p * points[i * 2 + 1];
        }
        return {x: x, y: y};
    }

    /**
     * 判断点是否在n阶贝塞尔曲线上
     * @param px 点的x坐标
     * @param py 点的y坐标
     * @param points 贝塞尔曲线的点集
     * @param tolerance 误差范围 最小为0.00001 实际能不能达到与迭代的次数有关
     * @param maxIteration 最大迭代次数 最小为1
     * @returns {Point | undefined} 点是否在三次贝塞尔曲线上
     */
    export function pointOnBezierCurve(px: number, py: number, points: number[], tolerance: number = 0.00001, maxIteration: number = 100): Point | undefined {
        // n 阶贝塞尔曲线的公式 B(t) = Σ C(n, i) * (1 - t)^(n - i) * t^i * P(i)
        // 其中
        // n 为 n 阶
        // i 为 0 到 n
        // C(n, i) 为 组合数 C(n, i) = n! / i! * (n - i)!
        // P(i) 为 控制点
        // t 为 0 到 1
        // B(t) 为 曲线上的点
        // 如果点不在多个点组成的矩形范围内 则直接返回
        let minX = points[0];
        let minY = points[1];
        let maxX = points[0];
        let maxY = points[1];
        let count = 0;
        for (let i = 2; i < points.length; i += 2, count++) {
            minX = Math.min(minX, points[i]);
            minY = Math.min(minY, points[i + 1]);
            maxX = Math.max(maxX, points[i]);
            maxY = Math.max(maxY, points[i + 1]);
        }
        if (px < minX || px > maxX || py < minY || py > maxY) {
            return undefined;
        }
        if (tolerance < 0.00001) {
            tolerance = 0.00001;
        }
        if (maxIteration < 1) {
            maxIteration = 1;
        }
        const enterPoint = {x: px, y: py};
        const step = 0.01;
        const distanceFunc = function (t: number) {
            const tp = bezierCurve(t, points, count);
            return Math.sqrt(Math.pow(enterPoint.x - tp.x, 2) + Math.pow(enterPoint.y - tp.y, 2));
        }
        let minDistance = Number.MAX_VALUE;
        let minDistanceT: number = 0;
        for (let t = 0; t <= 1; t += step) {
            const distance = distanceFunc(t);
            if (distance < minDistance) {
                minDistance = distance;
                minDistanceT = t;
            }
        }
        let minT = Math.min(0, minDistanceT - step);
        let maxT = Math.max(1, minDistanceT + step);
        // 使用二分法让距离减小到误差范围内
        for (let i = 0; i < maxIteration; i++) {
            const t = (minT + maxT) / 2;
            const distance = distanceFunc(t);
            if (distance < minDistance) {
                minDistance = distance;
                minDistanceT = t;
                minT = t;
            } else {
                maxT = t;
            }
            if (Math.abs(minDistance) <= tolerance) {
                break;
            }
        }
        return bezierCurve(minDistanceT, points, count);
    }
}
