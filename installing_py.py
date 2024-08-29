
def read_numbers_and_compute_hypotenuse_square(file_path: object) -> object:
    # Read the file containing two numbers
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Process the lines to extract numbers, assuming they might be on the same line or on separate lines
    if len(lines) == 1:
        a, b = map(int, lines[0].split())
    else:
        a = int(lines[0].strip())
        b = int(lines[1].strip())

    # Compute the square of the hypotenuse using the previously defined function
    return hypotenuse_square(a, b)

 read_numbers_and_compute_hypotenuse_square("/Users/zxj/Documents/Rosalindpractice/rosalind_ini2.txt")

