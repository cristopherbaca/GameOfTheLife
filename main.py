import numpy as np
import pygame

CELL_SIZE = 6 
SPEED = 10
FPS = 60
RANDOM_DENSITY = 0.25 
WHITE, BLACK = (255, 255, 255), (0, 0, 0)

def count_neighbors(board):
    neighbors = np.zeros(board.shape, dtype=int) 

    for dy in [-1, 0, 1]: 
        for dx in [-1, 0, 1]: 
            if dy != 0 or dx != 0: 
                neighbors += np.roll(np.roll(board, dy, axis=0), dx, axis=1)
    return neighbors

def next_generation(board):
    neighbors = count_neighbors(board) 
    rows, cols = board.shape 
    new_board = np.zeros(board.shape, dtype=int)

    for row in range(rows):
        for col in range(cols):
            if board[row, col] == 1: 
                if 2 <= neighbors[row, col] <= 3:
                    new_board[row, col] = 1
            else: 
                if neighbors[row, col] == 3:
                    new_board[row, col] = 1 
    return new_board

def random_board(rows, cols):
    new_board = np.zeros((rows, cols), dtype=int) 

    for r in range(rows):
        for c in range(cols): 
            new_board[r, c] = 1 if np.random.random() < RANDOM_DENSITY else 0 
    return new_board 

def compute_layout(window_size):
    width, height = window_size
    cols, rows = width // CELL_SIZE, height // CELL_SIZE 
    if cols < 1: 
        cols = 1
    if rows < 1:
        rows = 1

    offset_x = (width - (cols * CELL_SIZE)) // 2 
    offset_y = (height - (rows * CELL_SIZE)) // 2  
    return rows, cols, (offset_x, offset_y)

def resize_board(board, new_rows, new_cols):
    new_board = np.zeros((new_rows, new_cols), dtype=int) 
    for row in range(min(new_rows, board.shape[0])): 
        for col in range(min(new_cols, board.shape[1])):
            new_board[row, col] = board[row, col] 
    return new_board 

def create_window():
    desktop_h, desktop_w = pygame.display.get_desktop_sizes()[0]

    width = int(desktop_w * 0.9)
    height = int(desktop_h * 0.85) 

    pygame.display.set_mode((width, height), pygame.RESIZABLE)
    pygame.display.set_caption("Juego de la Vida")

    try:
        from pygame._sdl2 import Window
        window = Window.from_display_module() 
        #window = pygame.Window(title="Juego de la Vida", size=(width, height), resizable=True).maximize() 
        window.maximize() 
    except Exception:
        pass 

def draw_board(screen, board, offset):
    screen.fill(BLACK) 
    off_x, off_y = offset 
    rows, cols = board.shape 

    for row in range(rows): 
        for col in range(cols):
            if board[row, col] == 1: 
                x = off_x + (col * CELL_SIZE) 
                y = off_y + (row * CELL_SIZE) 
                rect = pygame.Rect(x, y, CELL_SIZE, CELL_SIZE) 
                pygame.draw.rect(screen, WHITE, rect) 

pygame.init() 
create_window() 
clock = pygame.time.Clock() 
window_size = pygame.display.get_surface().get_size() 
rows, cols, offset = compute_layout(window_size) 
board = np.zeros((rows, cols), dtype=int) 
paused = True 
generation = 0
elapsed_ms = 0
running = True

while running: 
    dt = clock.tick(FPS) 
    screen = pygame.display.get_surface()  
    if screen.get_size() != window_size: 
        window_size = screen.get_size() 
        rows, cols, offset = compute_layout(window_size) 
        board = resize_board(board, rows, cols) 
    for event in pygame.event.get(): 
        if event.type == pygame.QUIT: 
            running = False 
        elif event.type == pygame.KEYDOWN: 
            if event.key == pygame.K_SPACE: 
                paused = not paused 
                elapsed_ms = 0
            elif event.key == pygame.K_r and paused: 
                board = random_board(rows, cols) 
                generation = 0 
    if not paused: 
        elapsed_ms += dt 
        interval = 1000 / SPEED
        while elapsed_ms >= interval: 
            board = next_generation(board) 
            generation += 1 
            elapsed_ms -= interval 

    draw_board(screen, board, offset) 
    pygame.display.flip() 

pygame.quit()
