import numpy as np
import pygame

CELL_SIZE = 6
SPEED = 10
FPS = 60
RANDOM_DENSITY = 0.25 
WHITE, BLACK = (255, 255, 255), (0, 0, 0)

def count_neighbors(board):
    """
    Cuenta los vecinos vivos de cada celula usando un tablero toroidal para contarlos.
    """
    neighbors = np.zeros(board.shape, dtype=int) 

    for dy in [-1, 0, 1]: 
        for dx in [-1, 0, 1]: 
            if dy != 0 or dx != 0: 
                neighbors += np.roll(np.roll(board, dy, axis=0), dx, axis=1)
    return neighbors

def next_generation(board):
    neighbors = count_neighbors(board) 
    alive = board == 1
    survives = alive & ((neighbors == 2) | (neighbors == 3))
    born = ~alive & (neighbors == 3)

    return (survives | born).astype(int)

def random_board(rows, cols):
    return (np.random.random((rows, cols)) < RANDOM_DENSITY).astype(int)

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
    rows = min(new_rows, board.shape[0])
    cols = min(new_cols, board.shape[1])
    new_board[:rows, :cols] = board[:rows, :cols]
    return new_board 

def draw_point(event, board, offset, rows, cols):
    x, y = event.pos
    col = (x - offset[0]) // CELL_SIZE
    row = (y - offset[1]) // CELL_SIZE
    if 0 <= row < rows and 0 <= col < cols:
        board[row, col] = 1 - board[row, col]

def create_window():
    desktop_h, desktop_w = pygame.display.get_desktop_sizes()[0]
    width = int(desktop_w * 0.9)
    height = int(desktop_h * 0.85) 
    window = pygame.Window("Juego de la Vida", (width, height), resizable=True)

    # Se usa un try-except por si el manejador de ventanas no soporta el maximize().
    try:
        window.maximize()
    except pygame.error as error:
        print(f"Aviso: no se pudo maximizar la ventana ({error}).")

    return window

def draw_board(screen, board, offset):
    screen.fill(BLACK) 
    # Surfarray espera una matriz [x, y] y se tiene [y, x] por lo que se da la transpuesta de la matriz.
    gray = (board.T * 255).astype(np.uint8)
    rgb = np.repeat(gray[:, :, np.newaxis], 3, axis=2)
    # Se crea una superficie pequena de 1PX que luego se escala, eso mejora el rendimiento en vez de dibujar cada celula una por una.
    small = pygame.surfarray.make_surface(rgb)
    rows, cols = board.shape
    big = pygame.transform.scale(small, (cols * CELL_SIZE, rows * CELL_SIZE))
    screen.blit(big, offset)

class Game:

    def __init__(self):
        pygame.init()
        self.window = create_window()
        self.clock = pygame.time.Clock()
        self.window_size = self.window.get_surface().get_size()
        self.rows, self.cols, self.offset = compute_layout(self.window_size)
        self.board = np.zeros((self.rows, self.cols), dtype=int)
        self.paused = True
        self.generation = 0
        self.elapsed_ms = 0
        self.running = True

    def check_resize(self, screen):
        if screen.get_size() != self.window_size:
            self.window_size = screen.get_size()
            self.rows, self.cols, self.offset = compute_layout(self.window_size)
            self.board = resize_board(self.board, self.rows, self.cols)

    def handle_events(self):
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                self.running = False
            elif event.type == pygame.KEYDOWN:
                if event.key == pygame.K_SPACE:
                    self.paused = not self.paused
                    self.elapsed_ms = 0
                elif event.key == pygame.K_r and self.paused:
                    self.board = random_board(self.rows, self.cols)
                    self.generation = 0
            elif event.type == pygame.MOUSEBUTTONDOWN and event.button == 1 and self.paused:
                draw_point(event, self.board, self.offset, self.rows, self.cols)

    def update(self, dt):
        if self.paused:
            return

        # Si un frame tarda mucho, se calculan varias generaciones para no perder el ritmo.
        self.elapsed_ms += dt
        interval = 1000 / SPEED
        while self.elapsed_ms >= interval:
            self.board = next_generation(self.board)
            self.generation += 1
            self.elapsed_ms -= interval

    def draw(self, screen):
        draw_board(screen, self.board, self.offset)
        self.window.flip()

    def run(self):
        while self.running:
            dt = self.clock.tick(FPS)
            screen = self.window.get_surface()
            self.check_resize(screen)
            self.handle_events()
            self.update(dt)
            self.draw(screen)

        pygame.quit()

Game().run()
