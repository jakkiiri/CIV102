import math

TENSILE_STRENGTH, COMPRESSIVE_STRENGTH, SHEAR_STRENGTH_BOARD, POISSONS_RATIO, SHEAR_STRENGTH_GLUE, YOUNGS_MODULUS = 30, 6, 4, 0.2, 2, 4000
MAX_MOMENT, MAX_SHEAR = 69430, 257.3 # change this


def get_bounds(rectangles):
    if not rectangles:
        return -1, -1, 1, 1
    all_x = [coord for rect in rectangles for coord in [rect[0], rect[2]]]
    all_y = [coord for rect in rectangles for coord in [rect[1], rect[3]]]
    return min(all_x), min(all_y), max(all_x), max(all_y)

def calculate_centroid(rectangles):
    total_area, total_y = 0, 0
    for x1, y1, x2, y2 in rectangles:
        b, h = abs(x2 - x1), abs(y2 - y1)
        area, y_centroid = b * h, (y1 + y2) / 2
        total_area += area
        total_y += area * y_centroid
    return total_y / total_area if total_area else 0

def calculate_second_moment_of_area(rectangles, ybar):
    I_total = 0
    for x1, y1, x2, y2 in rectangles:
        b, h = abs(x2 - x1), abs(y2 - y1)
        area, y_centroid = b * h, (y1 + y2) / 2
        dy = y_centroid - ybar
        I_rect = (b * h**3) / 12
        I_total += I_rect + area * dy**2
    return I_total

def calculate_Q(rectangles, at_y, ybar):
    Q_total = 0
    for x1, y1, x2, y2 in rectangles:
        x_left = min(x1, x2)
        x_right = max(x1, x2)
        y_top = max(y1, y2)
        y_bottom = min(y1, y2)

        if y_top <= at_y:
            continue
        elif y_bottom >= at_y:
            y_clipped_bottom = y_bottom
        else:
            y_clipped_bottom = at_y

        h = y_top - y_clipped_bottom
        if h <= 0:
            continue

        b_rect = x_right - x_left
        area = b_rect * h
        y_centroid = (y_top + y_clipped_bottom) / 2
        dy = y_centroid - ybar

        Q_total += area * dy

    return Q_total

def calculate_b(rectangles, at_y):
    return sum(abs(x2 - x1) for x1, y1, x2, y2 in rectangles if y1 <= at_y <= y2 or y2 <= at_y <= y1)

def calculate_tau(Q_total, I, b_total):
    return MAX_SHEAR * Q_total / (I * b_total) if I and b_total else 0

def find_ybar_I_tau(rectangles):
    ybar = calculate_centroid(rectangles)
    I = calculate_second_moment_of_area(rectangles, ybar)
    Q_centroid = calculate_Q(rectangles, ybar, ybar)
    b_centroid = calculate_b(rectangles, ybar)
    Q_glue = calculate_Q(rectangles, 75, ybar)
    b_glue = calculate_b(rectangles, 75)
    tau_centroid = calculate_tau(Q_centroid, I, b_centroid)
    tau_glue = calculate_tau(Q_glue, I, 10)
    if I:
        min_x, min_y, max_x, max_y = get_bounds(rectangles)
        sigma_top = -MAX_MOMENT * (max_y - ybar) / I
        sigma_bottom = -MAX_MOMENT * (min_y - ybar) / I
    else:
        sigma_top = sigma_bottom = 0
    return ybar, I, Q_centroid, Q_glue, tau_centroid, tau_glue, sigma_bottom, sigma_top

def calculate_thin_plate_buckling(case, t, b):
    if case == 1:
        k = 4
    elif case == 2:
        k = 0.425
    else:
        k = 6
    return (k * (math.pi**2) * YOUNGS_MODULUS) * ((t / b)**2) / (12 * (1 - POISSONS_RATIO**2))

def calculate_thin_plate_buckling_shear(t, h, a):
    return (5 * (math.pi**2) * YOUNGS_MODULUS) * ((t / h)**2 + (t / a)**2) / (12 * (1 - POISSONS_RATIO**2))

def find_thin_plate_buckling(ybar):
    data = [[[0, 0, 10.635, 1.27, 2], [10.635, 0, 89.365, 1.27, 1], [89.365, 0, 100, 1.27, 2]], [[0, ybar, 1.27, 75, 3]], [[0, ybar, 1.27, 75, 3]], [], [], []]
    rect_count, region_count = -1, -1
    ans_rect, ans_region, ans_buckle = -1, -1, -1
    for rectangle in data:
        region_count = 0
        for region in rectangle:
            if region[4] == 3:
                current_buckle = calculate_thin_plate_buckling(3, region[2] - region[0], region[3] - region[1])
            else:
                current_buckle = calculate_thin_plate_buckling(region[4], region[3] - region[1], region[2] - region[0])
            if current_buckle > ans_buckle:
                ans_rect, ans_region, ans_buckle = rect_count, region_count, current_buckle
            region_count += 1
        rect_count += 1
    shear_buckling = calculate_thin_plate_buckling_shear(1.27, 75 - 1.27, 400)
    return ans_rect, ans_region, ans_buckle, shear_buckling
