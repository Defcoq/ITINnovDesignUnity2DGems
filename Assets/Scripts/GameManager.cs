using UnityEngine;
using UnityEngine.UI;
using UnityEngine.SceneManagement;
using System.Collections;

public class GameManager : MonoBehaviour {
  public static GameManager instance;

  private int score, highScore;
  public Text scoreText, highScoreText;

  private float time;
  public Text timeText;

  public bool started, gameOver;

  public GameObject gameOverPanel;
  public Text gameOverScore, gameOverHighScore;
  public Button playAgain;

  void Start() {
    instance = GetComponent<GameManager>();

    score = 0;
    scoreText.text = "Score: " + score;
    highScore = PlayerPrefs.GetInt("HighScore");
    highScoreText.text = "Score: " + highScore;

    time = 45;
    gameOver = false;
    UpdateTime();
  }

  void Update() {
    if (started) {
      time -= Time.deltaTime;
      UpdateTime();
      if (time <= 0) {
        GameOver();
      }
    }
  }

  public void IncreaseScore() {
    score++;
    scoreText.text = "Score: " + score;
    time += 5;
    UpdateTime();
  }

  public void StartTimer() {
    time = 45;
    started = true;
  }

  public void UpdateTime() {
    string minutes = Mathf.Floor(time / 60).ToString("00");
    string seconds = Mathf.Floor(time % 60).ToString("00");
    timeText.text = string.Format("Time: {0}:{1}", minutes, seconds);
  }

  private void GameOver() {
    time = 0;
    UpdateTime();
    started = false;
    gameOver = true;
    gameOverPanel.SetActive(true);
    gameOverScore.text = "Final Score: " + score;
    if (score > highScore) {
      gameOverHighScore.text = "High Score: " + score;
      PlayerPrefs.SetInt("HighScore", score);
    } else {
      gameOverHighScore.text = "High Score: " + score;
    }
  }

  public void LoadGame() {
    SceneManager.LoadScene("GameScene");
  }
}